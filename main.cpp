#include <iostream>
#include <fstream>
#include <set>
#include <array>
#include <algorithm>
#include <iomanip>
#include <mutex>
#include <ultimaille/all.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

template<typename T>
inline bool set_contains(const std::set<T>& set, T value) { return std::find(set.begin(), set.end(), value) != set.end(); }

const int facet_ordering[4][3] = {
    // <=> facet id to vertex ids
    // ordering used by libigl
    {0,1,2}, //facet 0
    {0,1,3}, //facet 1
    {1,2,3}, //facet 2
    {2,0,3}  //facet 3
};

void compute_tetra_adjacency_matrix(const UM::Tetrahedra& mesh, UM::CellAttribute<std::array<int,4>>& cell_adjacency);

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
    int nb_vertices = 0, nb_tetrahedra = 0, nb_walls = 0, nb_blocks = 0;

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
    UM::CellAttribute<int> block_id(mesh,-1);

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

    //should be the end of file, except for outputs of MC3D where the walls are also written
    ifs >> nb_walls;

    if(nb_walls > 0) {
        // case of a .hexex file from MC3D https://github.com/HendrikBrueckler/MC3D
        std::cout << "/!\\ The input file contains walls.\n"
                  << "    Because Graphite cannot display cell facet attributes,\n"
                  << "    the program will group cells by block\n"
                  << "    and write block IDs instead." << std::endl;

        std::set<std::set<int>> wall_triangles;
        for(int w = 0; w < nb_walls; w++) {
            int vertex_1 = 0, vertex_2 = 0, vertex_3 = 0;
            double distance_to_origin = 0.0;
            ifs >> vertex_1 >> vertex_2 >> vertex_3 >> distance_to_origin;
            std::set<int> vertex_ids{vertex_1,
                                     vertex_2,
                                     vertex_3};//the set will sort them
            wall_triangles.emplace(vertex_ids);
            //should be the end of line
        }
        //should be the end of file

        //construct the tetra adjacency matrix
        UM::CellAttribute<std::array<int,4>> cell_adjacency(mesh,{-1,-1,-1,-1});
        compute_tetra_adjacency_matrix(mesh,cell_adjacency);

        UM::DisjointSet cells_grouping(nb_tetrahedra);//will group cells by block
        //TODO better parallelization
        std::mutex disjont_set_mutex;
        for(int t = 0; t < nb_tetrahedra; t++) { //for each cell

            // idea : if the triangle between the current cell i and its j^th neighbor
            // is not in the set of wall triangles, then i and its j^th neighbor are in the same block
            // -> merge them in the disjont set data structure
            // if the "id" of its j^th neighbor is -1, this is a boundary facet

            #ifdef USE_OPENMP
            #pragma omp parallel for num_threads(4)//check the 4 facets in //
            #endif
            for(int f=0; f<4; f++) { //for each facet of the current cell
                if(cell_adjacency[t][f]!=-1) {
                    std::set<int> vertex_ids{mesh.cells[4*t+facet_ordering[f][0]],
                                             mesh.cells[4*t+facet_ordering[f][1]],
                                             mesh.cells[4*t+facet_ordering[f][2]]};//the set will sort them
                    if(!set_contains<std::set<int>>(wall_triangles,vertex_ids)) {//if this triangle is not a wall
                        std::lock_guard<std::mutex> guard(disjont_set_mutex);
                        cells_grouping.merge(t,cell_adjacency[t][f]);//the 2 cells are in the same block
                    }
                }
                //else it is a boundary facet
            }

            if( t%500==0 ) { //every 500 cells
                std::cout << "\rComputation of block IDs... " << std::setprecision(3) << (((double)t)/nb_tetrahedra)*100 << "%   " << std::flush;
            }
        }
        std::cout << std::endl;

        //get block ids
        nb_blocks = cells_grouping.get_sets_id(block_id.ptr.get()->data);

        UM::write_geogram(
            output_filename,    //file name
            mesh,               //volumetric mesh
            UM::VolumeAttributes{
                {},                                         //per point attributes
                { {"block_id", block_id.ptr} },             //per cell attributes
                {},                                         //per cell facet attributes
                { {"u", u.ptr}, {"v",v.ptr}, {"w",w.ptr} }  //per cell corner attributes
            }
        );
    }
    else {
        //write without block_id
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
    }

    std::cout << "Done. They are:\n"
              << nb_vertices << " vertices\n"
              << nb_tetrahedra << " tetrahedra" << std::endl;
    if(nb_walls > 0) {
        std::cout << nb_walls << " wall triangles\n"
                  << nb_blocks << " blocks" << std::endl;
    }

    return 0;
}

// WARNING cell_adjacency must be filled with -1
void compute_tetra_adjacency_matrix(const UM::Tetrahedra& mesh, UM::CellAttribute<std::array<int,4>>& cell_adjacency) {

    // inspired by tet_tet_adjacency() of libigl
    // https://libigl.github.io/
    // https://github.com/libigl/libigl/blob/main/include/igl/tet_tet_adjacency.cpp
    // this part of the code is under the Mozilla Public License v2.0 http://mozilla.org/MPL/2.0/

    std::vector<std::array<int,5>> TTT(4*mesh.ncells());//(4*nb_tetrahedra = nb cell facets) lines, 5 columns

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for(int t=0; t<mesh.ncells(); t++) { //for each tetra
        for(int f=0; f<4; f++) { //for each facet of the current tetra
            TTT[4*t+f] = {
                //store 5 elements : the 3 vertices, the cell index and the facet index
                mesh.cells[4*t+facet_ordering[f][0]],//first vertex
                mesh.cells[4*t+facet_ordering[f][1]],//second vertex
                mesh.cells[4*t+facet_ordering[f][2]],//third vertex
                t,
                f
            };
            //sort vertex ids
            std::sort(TTT[4*t+f].begin(), TTT[4*t+f].begin()+3);//the 3 first elements
        }
    }
    
    std::sort(TTT.begin(),TTT.end());//sort matrix by first column = by lowest vertex id
    
    for(int i=1; i<TTT.size(); i++) {//for each cell facet
        const std::array<int,5>& r1 = TTT[i-1];//previous facet
        const std::array<int,5>& r2 = TTT[i];  //current facet
        if((r1[0]==r2[0]) && (r1[1]==r2[1]) && (r1[2]==r2[2])) {//if the 2 facets have the same 3 vertices
            cell_adjacency[r1[3]][r1[4]] = r2[3];//link cell of previous facet to cell of current facet
            cell_adjacency[r2[3]][r2[4]] = r1[3];//link cell of current facet to cell of previous facet
        }
    }
}