#include <iostream>
#include <fstream>
#include <ultimaille/all.h>

int main(int argc, char *argv[]) {

    if(argc < 3) {
        std::cerr << "Wrong usage, it should be:" << std::endl;
        std::cerr << "\t hexex2geogram <input.hexex> <output.geogram>" << std::endl;
        return 1;
    }

    std::string input_filename  = argv[1],
                output_filename = argv[2];

    std::ifstream ifs(input_filename);
    if(!ifs.is_open()) {
        std::cerr << "Unable to open '" << input_filename << "'" << std::endl;
        return 1;
    }

    UM::Tetrahedra mesh;
    int nb_vertices = 0, nb_tetrahedra = 0;

    // read number of vertices
    ifs >> nb_vertices;
    mesh.points.resize(nb_vertices);

    // read vertices
    for(int n = 0; n < nb_vertices; n++) {
        double x,y,z;
        ifs >> x >> y >> z;
        UM::vec3 coordinates(x,y,z);
        mesh.points[n] = coordinates;
    }

    // read number of tetrahedra
    ifs >> nb_tetrahedra;
    mesh.cells.resize(4*nb_tetrahedra);
    UM::CellCornerAttribute<double> u(mesh,nan("")),
                                    v(mesh,nan("")),
                                    w(mesh,nan(""));

    // read cells (tetrahedra)
    for(int t = 0; t < nb_tetrahedra; t++) {
        //read vertex indices of cell corners
        for(int j = 0; j < 4; j++) {
            int vertex_index = 0;
            ifs >> vertex_index;
            mesh.cells[4*t+j] = vertex_index;
        }
        //read param at corner 0, 1, 2, 3
        for(int v_id = 0; v_id < 4; v_id++) {
            //read u, v, w values for this cell corner
            double current_u, current_v, current_w;
            ifs >> u[4*t+v_id] >> v[4*t+v_id] >> w[4*t+v_id];
        }
        //should be the end of line
    }
    //should be the end of file, except for outputs of MC3D where the walls are alsy written

    std::cout << "there are " << mesh.nverts() << " vertices and " << mesh.ncells() << " tetrahedra" << std::endl;
    
    UM::write_geogram(
        output_filename,    //file name
        mesh,               //volumetric mesh
        UM::VolumeAttributes{
            {},                                         //per point attributes
            {},                                         //per cell attributes
            {},                                         //per cell facet attributes
            { {"u", u.ptr}, {"v",v.ptr}, {"w",w.ptr} }  //per cell corner attributes
        }
    );

    return 0;
}