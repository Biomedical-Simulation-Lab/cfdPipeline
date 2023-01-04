""" Generate volume mesh and relevant submission files.
"""

from pathlib import Path 
import pyvista as pv 
from geometry_tools.meshing import Mesher
from geometry_tools.make_submission_file import SubmissionTemplate
import sys
import time
from datetime import timedelta

def volume_meshing(proj_dir):
    proj_dir = Path(proj_dir)
    clipped_dir = proj_dir / '02_processed' 
    surf_file = sorted(clipped_dir.glob('*.vtp'))[0]

    proj_dir = surf_file.parents[1]

    points_dir = proj_dir / '02_points' 

    point_file = points_dir / (surf_file.stem + '_endpoints.vtm') 
    assert point_file.exists()

    mesh_out_dir = proj_dir / '03_mesh' 
    data_out_dir = proj_dir / '03_data'  
    submission_out_dir = proj_dir / '03_submissions' 

    for f in [mesh_out_dir, data_out_dir, submission_out_dir]:
        if not f.exists():
            f.mkdir(parents=True, exist_ok=True)

    vtufile = mesh_out_dir / (surf_file.stem + '.vtu')
    meshfile = data_out_dir / (surf_file.stem + '.h5')
    infofile = data_out_dir / (surf_file.stem + '.info')
    xmlgzfile = data_out_dir / (surf_file.stem + '.xml.gz')

    start = time.time()
    print('\n' + surf_file.stem)

    if not xmlgzfile.exists():
        # try:
        surf = pv.read(surf_file)

        points = pv.read(point_file)
        in_points = points['inlets'].points
        out_points = points['outlets'].points
        an_points = points['aneurysms'].points

        m = Mesher(
            surf,
            inlet_points=in_points, 
            outlet_points=out_points,
            aneurysm_points=an_points,
            )
        m.update_inlets_outlets()

        if not vtufile.exists():
            m.generate_volume_mesh()
            m.mesh.save(vtufile)
        else:
            m.mesh = pv.read(vtufile)

        m.update_inlets_outlets()

        m.generate_centerlines()
        m.generate_centerlines(include_aneurysms=False)

        m.generate_flow_rates()
        m.generate_h5_file(meshfile)
        m.generate_flow_rates_legacy()
        m.update_inlets_outlets()
        m.generate_info_file(infofile,)  
        m.generate_xml_gz_file(xmlgzfile)

        # Create submission file
        s = SubmissionTemplate(meshfile.stem)
        s.save_script(submission_out_dir)

        print('\n Case done', timedelta(seconds=time.time() - start))

    else:
        print('Output file exists.')



if __name__ == "__main__":
    proj_dir = Path(sys.argv[1]) 
    volume_meshing(proj_dir=proj_dir)
