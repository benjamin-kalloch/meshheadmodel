#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <CGAL/ImageIO/analyze.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/odt_optimize_mesh_3.h>
#include <CGAL/perturb_mesh_3.h>

#include <vtkUnstructuredGrid.h>    
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>       // istringstream
#include <chrono>        // for processing time measurements
#include <cstdio>        // std::remove
#include <getopt.h>      // use GNU extension of getopt to be able to parse more readable multi-character options
#include <vector>

#include "Gmsh.h"
#include "GModel.h"
#include "GEntity.h"
#include "MTetrahedron.h"
#include "MTriangle.h"
#include "MElement.h"


typedef int Curve_index;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    // labeled multi-domain 3D image:
    // This class includes a member function that provides, by interpolation, the subdomain index of any query point.
    // An intersection between a segment and bounding surfaces is detected when both segment endpoints are associated
    // with different values of subdomain indices. The intersection is then constructed by bisection. The bisection stops
    // when the query segment is shorter than a given error bound e. This error bound is given by e=d ×bound where d is
    // the length of the diagonal of the bounding box (in world coordinates) and bound is the argument passed to the
    // constructor of Labeled_image_mesh_domain_3.
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,Kernel> Image_domain;

    // traits class providing the types and functors required to implement the intersection tests and intersection
    // computations for polyhedral boundary surfaces.
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> Polyhedron_domain;
typedef CGAL::Mesh_polyhedron_3<Kernel>::type Polyhedron;

    // To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// function prototypes
void printHelpText();
bool createPolyhedronFromOFF( const char*, Polyhedron& );

//
// Helper class that allows the use of a compound domain 
// consisting of one Labeled_image_domain and multiple Polyhedral_domain(s).
// 
class Hybrid_domain
{
  const Image_domain& m_image_domain;                   // domain containing the segmentation image
  const std::vector<Polyhedron_domain*> m_poly_domains; // vector of domains for an unknown number of electrode surfaces

public:
  Hybrid_domain( const Image_domain& _image_domain, const std::vector<Polyhedron_domain*> _poly_domains) :
    m_image_domain( _image_domain ),
    m_poly_domains( _poly_domains )
  {}

  // types required by the 'MeshDomain_3' concept
  typedef int Surface_patch_index;
  typedef int Subdomain_index;
  typedef int Index;

  typedef Kernel R;
  typedef Kernel::Point_3 Point_3;
  typedef CGAL::cpp11::tuple<Point_3, Index, int> Intersection;

  //
  // The bounding box of the Hybrid_domain is the union of bounding-boxes of the subdomains.
  //
  CGAL::Bbox_3 bbox() const {
    CGAL::Bbox_3 hybrid_domain_bbox( m_image_domain.bbox() );

    // add electrode Bboxes to Label_image_domain Bbox
    std::for_each( m_poly_domains.begin(),
                   m_poly_domains.end(),
                   [&hybrid_domain_bbox] (Polyhedron_domain* _domain) { 
                        hybrid_domain_bbox + _domain->bbox();
                 });
    
    return hybrid_domain_bbox;
    }

  //
  // Initiate seeding in all subdomains.
  //
  struct Construct_initial_points
  {
    Construct_initial_points(const Hybrid_domain& _domain) :
        r_domain_( _domain ) {}
    
    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 100) const
    {
        int numPointsToConstructPerSubdomain = n / ( r_domain_.m_poly_domains.size() + 1 );
        // due to integer rounding n != totalPointToConstruct might be possible
        int totalPointToConstruct = numPointsToConstructPerSubdomain * ( r_domain_.m_poly_domains.size() + 1 );

        // construct initial seed points on image domain
        typedef Image_domain::Index Image_Domain_Index;
        std::vector< std::pair<Point_3, Image_Domain_Index> > image_domain_points_vector;
        Image_domain::Construct_initial_points cstr_image_domain_initial_points 
            = r_domain_.m_image_domain.construct_initial_points_object();  // assign the remaining points to the image domain
        cstr_image_domain_initial_points( std::back_inserter( image_domain_points_vector ), numPointsToConstructPerSubdomain + (n - totalPointToConstruct));

        for( std::size_t i = 0, end = image_domain_points_vector.size(); i < end; i++ ) 
        {
            *pts++ = std::make_pair( image_domain_points_vector[i].first, 2 );
        }
        // construct inital seed points on all polyhedral domains
        for( std::size_t d = 0; d < r_domain_.m_poly_domains.size(); d++)
        {
            Polyhedron_domain *poly_domain = r_domain_.m_poly_domains.at(d);

            typedef Polyhedron_domain::Index Poly_Domain_Index;
            std::vector< std::pair< Point_3, Poly_Domain_Index > > poly_domain_points_vector;
            Polyhedron_domain::Construct_initial_points cstr_poly_domain_initital_points
                = poly_domain->construct_initial_points_object();    
            cstr_poly_domain_initital_points( std::back_inserter( poly_domain_points_vector ), numPointsToConstructPerSubdomain );
    
            for( std::size_t i = 0, end = poly_domain_points_vector.size(); i < end; i++ )
            {
                *pts++ = std::make_pair( poly_domain_points_vector[i].first, 1 );
            }
        }

        return pts;
    }
    private:
        const Hybrid_domain& r_domain_;
  }; // end Construct_initial_points

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points( *this );
  }

  //
  // Check inclusion if query point for all subdomains. 
  //
  struct Is_in_domain
  {
    Is_in_domain( const Hybrid_domain& domain ) :
        r_domain_( domain ) 
    {}

    boost::optional< Subdomain_index > operator()(const Kernel::Point_3& p) const
    {
        boost::optional<Subdomain_index> subdomain_index = 
            r_domain_.m_image_domain.is_in_domain_object()(p);
        
        // contains a value if we are in a subdomain of the image domain
        if( subdomain_index ) 
        {
            return subdomain_index;
        }
        // if not in image domain check if in polyhedral_domain
        else
        {
            for(std::size_t d = 0; d < r_domain_.m_poly_domains.size(); d++)
            {
                Polyhedron_domain *poly_domain = r_domain_.m_poly_domains.at(d);

                boost::optional<Subdomain_index> poly_domain_index = 
                    poly_domain->is_in_domain_object()(p);

                // return index of polyhedral domain (+100 to separate indices from image-subdomain indices)
                if( poly_domain_index )
                {
                    return 100 + d;
                }
            }
        }

        // return empty subdomain_index indicating that the query point is neither in the
        // Label_image_domain nor in one of the Polyhedron_domains
        return subdomain_index;
    }

    private:
    const Hybrid_domain& r_domain_;
  }; // end Is_in_domain
  
  Is_in_domain is_in_domain_object() const { return Is_in_domain( *this ); }

  //
  // Check for intersection with the domain boundaries.
  // For the Image_domain the voxel coordinates are redturned (?).
  // For the Polyhedron_domain an interpolated coordinate (over the surface) is returned.
  //
  struct Construct_intersection
  {
    Construct_intersection(const Hybrid_domain& domain) :
        r_domain_(domain)
    {}

    template <typename Query>    // Query can be a Segment_3, Ray_3, Line_3
    Intersection operator()(const Query& query) const
    {
        using boost::get;
        
        // intersection with image domain
        Image_domain::Intersection image_domain_intersec =
            r_domain_.m_image_domain.construct_intersection_object()(query);

        // if intersection was found, return it 
        if( get<2>(image_domain_intersec) != 0)
        {
            return Intersection(get<0>(image_domain_intersec), 2, get<2>(image_domain_intersec));
        }


        for(std::size_t d = 0; d < r_domain_.m_poly_domains.size(); d++)
        {
            Polyhedron_domain *poly_domain = r_domain_.m_poly_domains.at(d);

            // interpolation with polyhedral domain
            Polyhedron_domain::Intersection poly_domain_intersec =
                poly_domain->construct_intersection_object()(query);

            // if intersection was found return it
            if( get<2>(poly_domain_intersec) != 0)
            {
                const Point_3 inter_point = get<0>(poly_domain_intersec);
                
                if( !r_domain_.m_image_domain.is_in_domain_object()(inter_point) )
                {    
                    return Intersection( inter_point, 1, get<2>(poly_domain_intersec) );
                }
            }
        }

        // no intersection found
        return Intersection();
    }
      private:
        const Hybrid_domain& r_domain_;
  }; // end Construct_intersection

  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection( *this );
  }

  // Index types converters
  Index index_from_surface_patch_index( const Surface_patch_index& index) const
  {
    return index;
  }

  Index index_from_subdomain_index( const Subdomain_index& index) const
  {
    return index;
  }

  Surface_patch_index surface_patch_index(const Index& index) const
  {
    return index;
  }    

  Subdomain_index subdomain_index(const Index& index) const
  {
    return index;
  }
};    // end class Hybrid_domain

//
// Wrapper domain for the Hybrid_domain.
// The PolylineFeatures allow sharp edges at the electrode corners,
// however the location of these edges must be specified in a
// seperate text file.
// 
typedef CGAL::Mesh_domain_with_polyline_features_3<Hybrid_domain> Domain;

//
// Define the type of triangulation
//
typedef CGAL::Mesh_triangulation_3<Domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

//
// Mesh criteria
//
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;


int main(int argc, char*argv[])
{
  // ******* parsing of commandline arguments **********
  const char *imagefile      = nullptr,
             *outfile        = nullptr,
             *smoothskinfile = nullptr;

  std::vector<std::string> electrodepaths;
 
  float faceAngBound   = 30.f,
        faceSizeBound  = 6.f,
        faceDistBound  = 4.f,
        cellREratio    = 3.f,
        cellSizeBound  = 8.f;

  int use_lloyd = 0, use_odt = 0, perturb = 0, exude = 0;
  
  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"imagefile", 1, NULL, 0},        // name, has_arg (1=required, 2=optional), flag, val
    {"smoothskinfile", 1, NULL, 0},
    {"electrodefiles", 1, NULL, 0},
    {"outfile", 1, NULL, 0},
    {"f_angbound", 1, NULL, 0},
    {"f_sizebound", 1, NULL, 0},
    {"f_distbound", 1, NULL, 0},
    {"c_reratio", 1, NULL, 0},
    {"c_sizebound", 1, NULL, 0},
    {"help", 0, NULL, 0},
    {"lloyd",0,&use_lloyd, 1},
    {"odt",0,&use_odt, 1},
    {"perturb",0,&perturb, 1},
    {"exude",0,&exude, 1},
    {NULL, 0, NULL, 0}            // "The last element of the array must be filled with zeroes"
  };
                                            // no short options
  while( ( c = getopt_long_only( argc, argv, "", long_options, &option_index ) ) != -1 )
  {
    if( c == 0 )    // long option detected, c > 0 = short option detected, c = -1 end of options
    {
        switch( option_index )
        {
            case 0:    // imagefiles
                imagefile = optarg;    
            break;

            case 1: // smooth skin surface file
                smoothskinfile = optarg;
            break;

            case 2: // electrode files
            {
                std::istringstream electrodes_string( optarg );
                std::string path;
                for( std::size_t p = 0; std::getline(electrodes_string, path, ',' ); p++ )
                {
                    electrodepaths.push_back( path );
                }
            }
            break;

            case 3:    // outfile
                outfile = optarg;
            break;

            case 4:    // face ang-bound
                faceAngBound = std::stof(optarg);
            break;

            case 5:    // face size-bound
                faceSizeBound = std::stof(optarg);
            break;

            case 6:    // face distance-bound
                faceDistBound = std::stof(optarg);
            break;

            case 7:    // cell radius-edge ratio 
                cellREratio = std::stof(optarg);
            break;

            case 8:    // cell size-ratio
                cellSizeBound = std::stof(optarg);
            break;
            
            case 9:    // help
                printHelpText();
                 return 0;
            break;
        }
     }
   }

  if( !imagefile || !outfile ) 
  {
    std::cerr << "Imagefile and outfile must be specified!" << std::endl;
    return 1;
  }

 std::cout << "Starting volume mesh using the following parameters: \n"    
           << "\t imagefile = " << imagefile << "\n"
           << "\t smoothskinfile = " << ( smoothskinfile ? smoothskinfile : "none provided" ) << "\n"
           << "\t electrodes = ";
 std::for_each(electrodepaths.begin(), electrodepaths.end(), [](std::string &path){std::cout << path << " ";});
 std::cout << "\n\t outfile = " << outfile << "\n"
           << "\t faceAngBound = " << faceAngBound << "\n"
           << "\t faceSizeBound = " << faceSizeBound << "\n"
           << "\t faceDistBound = " << faceDistBound << "\n"
           << "\t cellREratio = " << cellREratio << "\n"
           << "\t cellSizeBound = " << cellSizeBound << "\n"
           << "\t use lloy-optimization = " << (use_lloyd == 1 ? "yes" : "no") << "\n"
           << "\t use odt-optimization = " << (use_odt == 1 ? "yes" : "no") << "\n"
           << "\t use perturbate mesh = " << (perturb == 1 ? "yes" : "no") << "\n"
           << "\t use exude mesh = " << (exude == 1 ? "yes" : "no") << std::endl;

  // ************* Read the input files *****************
    // 1) Load the segmentation image
    // possible input format:
    //     ANALYZE-image format (hdr + img file) in unsigned byte format, NOT short or unsigned int
    //  Why Byte? Because:
    //  By default, the class template `CGAL::Labeled_image_mesh_domain_3` only deals with
    //  labelized image whose data type is 8 bits (signed or not, that makes no difference). 
  CGAL::Image_3 image;
  if( !image.read(imagefile) )
  {
    std::cerr << "Error: Cannot read file " <<  imagefile << std::endl;
    return EXIT_FAILURE;
  }

    // 2) Load surfaces (electrodes and optionally smoothskin)
    // possible input format:
    //    OFF (object file format), both binary or ASCII
  size_t num_surface_files = electrodepaths.size() + ( smoothskinfile != nullptr );
  std::vector<Polyhedron_domain*> polydomains; // Polyherdron_domain as pointer to avoid copying
  std::vector<Polyhedron> polyhedrons( num_surface_files );  // we may have to represent skin as a polyhedron too, 
                                                            // but its features should not be preserved,
                                                            // i.e. it does not have any since its supposed to be *smooth*
  std::vector<Polyhedron>::iterator current_poly_it = polyhedrons.begin();
  
  // prepare data structures for feature detection
    // a vector of lists of corners of a poly_domain (one list per poly_domain)
  std::vector< std::vector< std::pair<Polyhedron_domain::Corner_index,Kernel::Point_3> > > corners( electrodepaths.size() );    
    // a vector of lists of curves of a poly_domain (one list per poly_domain)
  std::vector< std::vector< CGAL::cpp11::tuple<Curve_index, std::pair<Kernel::Point_3,Domain::Index>, std::pair<Kernel::Point_3,Domain::Index> > > > curves( electrodepaths.size() ); 

  // load smooth skin surface (if provided) 
  if( smoothskinfile )
  {
    if( ! createPolyhedronFromOFF( smoothskinfile, *current_poly_it ) ) return EXIT_FAILURE;
    polydomains.push_back( new Polyhedron_domain( *current_poly_it ) );
    current_poly_it++;
  }

  // load electrode surfaces
  for( std::size_t p = 0; current_poly_it != polyhedrons.end(); current_poly_it++, p++ )
  {
    if( ! createPolyhedronFromOFF( electrodepaths.at( p ).c_str(), *current_poly_it ) ) return EXIT_FAILURE;
    
    // create new Polyhedron_domain and store     // Polyhedron_domain has a call by reference construtcor
    
    polydomains.push_back( new Polyhedron_domain( *current_poly_it ) );
    // Detect sharp features and boundaries of the polyhedral components of the complex, and insert them as features of the domain.
    // Parameter: The maximum angle (in degrees) between the two normal vectors of adjacent triangles.
    //        For an edge of the polyhedron, if the angle between the two normal vectors of its
    //        incident facets is bigger than the given bound, then the edge is considered as a
    //        feature edge, and inserted as a feature of the domain. Defaults to 60.  
      polydomains.back()->detect_features( );
    // retrieve detected corners
    polydomains.back()->get_corners( std::back_inserter( corners.at( p ) ) );
    polydomains.back()->get_curves( std::back_inserter( curves.at( p ) ) );
  }


  // Image Domain:
  // A 3D labeled image is a grid of voxels, where each voxel is associated with an index (a subdomain index)
  // characterizing the subdomain in which the voxel lies. This class is a model of the concept MeshDomain_3.
  // The domain to be discretized is the union of voxels that have an non-default index (different from the
  // default constructed value of the type Image::Type).
  // This class includes a member function that provides, by interpolation, the subdomain index of any query
  // point. An intersection between a segment and bounding surfaces is detected when both segment endpoints
  // are associated with different values of subdomain indices. The intersection is then constructed by
  // bisection. The bisection stops when the query segment is shorter than a given error bound e. This error
  // bound is given by e=d*bound where 'd' is the length of the diagonal of the bounding box (in world 
  // coordinates) and 'bound' is the argument passed to the constructor of Labeled_image_mesh_domain_3.
  // Second parameter:
  // The default value for that second parameter (i.e. the rror bound) is 1e-3, and that is too low for
  // precise images like the one you use (256x256x256). That parameter sets the precision of the
  // bissection algorithm used to compute intersections with the implicit surfaces
  // of the image and query segments.
  // 'e' is relative to the size of the domain:
  // "This error bound is given by e=d×bound where d is the length of the diagonal of the bounding box"
  //  Its (absolute) value should stay smaller than facet_size or cell_size and facet_distance. 
  Image_domain image_domain(image, 1e-10);

  // the Hybrid_domain
  Domain domain( image_domain, polydomains );
  std::vector< std::vector< Kernel::Point_3 > > featured_curves;
  Curve_index last_curve_index;
  size_t p;
  std::vector<Polyhedron_domain*>::iterator polydomain_it;
  for
  ( p = 0, polydomain_it = ( smoothskinfile ? ++(polydomains.begin()) : polydomains.begin() ); // skip smooth-skin poly domain
    polydomain_it != polydomains.end();
    p++, polydomain_it++
  )
  {    
    // a vector of lists per domain of all of its curves
    std::vector< CGAL::cpp11::tuple<Curve_index, std::pair<Kernel::Point_3,Domain::Index>, std::pair<Kernel::Point_3,Domain::Index> > > &curves_of_domain = curves.at( p );
    
    // loop over all detected curves
    for( size_t l = 0 ; l < curves_of_domain.size(); l++)
      {
        std::vector< Kernel::Point_3> curve_polyline;    // a curve is a vector of Point_3s of the same curve index
        CGAL::cpp11::tuple<Curve_index, std::pair<Kernel::Point_3,Domain::Index>, std::pair<Kernel::Point_3,Domain::Index> > &current_curve = curves_of_domain.at( l );
        Curve_index ci = std::get<0>( current_curve );    // each element returned by 'get_curves' is a triple: (curve index, start pt of curve, end point of curve)
                                                        // start & end point are pairs of Point3 and the associated point index
                                                        // (i.e. to access Point3 object of the start point: std::get<1>( curves_of_domain.at( l ) ).first 
        
        // create points on the detected curve (equidistant 1 unit of measurement, but at least one point between start & end)
        Kernel::Point_3 startPoint = std::get<1>( current_curve ).first;
        Kernel::Point_3 endPoint   = std::get<2>( current_curve ).first;
        float curve_length = (*polydomain_it)->curve_length( ci );
        float step_size    = std::min(curve_length * 0.5f, 1.f);
        //std::cout << "Electrode #" << p << ", curve length = " << curve_length << std::endl;

        curve_polyline.push_back( startPoint );    // since I found no way of obtainig the vertices of the curve, we sample the curve between the start and end point
        for( float s = step_size; s < curve_length; s += step_size )
        {
            curve_polyline.push_back((*polydomain_it)->construct_point_on_curve( startPoint, ci, s ) );
        }
        curve_polyline.push_back( endPoint );

        featured_curves.push_back( curve_polyline );
    }
  }
  domain.add_features( featured_curves.begin(), featured_curves.end() );

  // ************** Meshing Generation *******************
  std::chrono::high_resolution_clock::time_point t_start; 
  std::chrono::high_resolution_clock::time_point t_end;
  auto tic = []()->std::chrono::high_resolution_clock::time_point { return std::chrono::high_resolution_clock::now(); };
  auto printDuration = [](std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)->void {
        std::cout << "Finished after " << std::chrono::duration_cast<std::chrono::seconds>( end - start ).count() << " seconds." << std::endl;
  };

  std::cout << "Creating raw volume mesh... " << std::flush;
  // Volume meshing according to provided mesh criteria:
  //    1) Note that each size or distance parameter can be specified using two ways: either as a scalar
  //       field or as a numerical value when the field is uniform.
  //       With the scalar field you can set the criteria space varying, with the scalar it will be
  //       uniform throughout the volume.
  //    2) Each parameter has a special default value 'ignored' which means that the corresponding criterion
  //       will be ignored.
  //    3) Additional parameters not used here:
  //        I)  edge_size : [field/scalar] upper bound for the lengths of curve segment edges.
  //                         This parameter has to be set to a positive value when 1-dimensional
  //                          features protection is used.
  //        II) facet_topology : parameter of type 'Mesh_facet_topology';  the set of topological
  //                             constraints which have to be verified by each surface facet. The
  //                             default value is CGAL::FACET_VERTICES_ON_SURFACE. 
  //                             See Mesh_facet_topology manual page to get all possible values
  // faceAngBound    ... a lower bound for the angles (in degrees) of the surface mesh facets            
  //                    Defines the miminum angle of a surface triangle. For angular_bound, a
  //                    value of 60 would mean a uniformly sized cell. The meshing algorithm is
  //                    proved to terminate if the lower bound on facets angles is not bigger
  //                    than 30 degrees
  // faceSizeBound    ... [field/scalar] an upper-bound or for the radii of the surface Delaunay balls
  // faceDistBound  ... [field/scalar] an upper bound for the absolute distance between the facet
  //                    circumcenter and the center of its surface Delaunay ball.
  //                    Controls the approximation error from mesh surface facets to the surface of
  //                    the domain you are meshing. A good compromise between acurate surface aproximation
  //                    and smoothness is sqrt(2) * voxelsize. If this value is too small the mesh will
  //                    follow the voxelized structure of the mesh, if too large the surfaces will be only
  //                    poorly approximated.
  // cellREratio    ... an upper bound for the radius-edge ratio of the mesh tetrahedra
  //                    This parameter controls the shape of mesh cells. The radius edge ratio of a
  //                    simplex (triangle or tetrahedron) is the the ratio between its circumradius
  //                     and its shortest edge. There is a theoretical bound for this parameter:
  //                the Delaunay refinement process is guaranteed to terminate for values of
  //                    radius_edge_ratio bigger than 2.
  //                    -> The smaller the value the better the shape of the cell
  // cellSizeBound  ... [field/scalar] an upper-bound for the circumradii of the mesh tetrahedra
  Mesh_criteria criteria(edge_size=faceSizeBound,   // previously set to 0.5
                         facet_angle=faceAngBound,
                         facet_size=faceSizeBound,
                         facet_distance=faceDistBound,
                         cell_radius_edge_ratio=cellREratio,
                         cell_size=cellSizeBound);    
  t_start = tic();
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, features(), no_perturb(), no_exude(), no_lloyd(), no_odt());
  t_end = tic();
  printDuration( t_start, t_end );

  // * Optimization phase * 
  // 1) the Lloyd and odt-smoother:
  // =  are global optimizers moving the mesh vertices to minimize a mesh energy
  // In both cases the mesh energy is the L1 error resulting from the interpolation of the function
  // f(x)=x^2 by a piecewise linear function.
  // The Lloyd optimizer is known to be blind to the occurrence of slivers in the mesh while the
  // odt-smoother tends to chase them out.
  // Both of them are global optimizers, meaning that they try to improve the whole mesh rather than
  // focusing on the worst elements. However, both are empirically known to be very efficient as a
  // preliminary step of optimization, as they tend to enhance the efficiency of the perturber and/or
  // exuder applied.
  // 2) the pertuber and exuder:
  // = focus on improving the worst mesh elements
  // The perturber improves the meshes by local changes in the vertices positions aiming to
  // make sliver disappear. 
  // The exuder chases the remaining slivers by re-weighting mesh vertices with optimal weights.
  //
  // Each optimization process can be activated or not, and tuned according to the user requirements
  // and the available time. By default, only the perturber and the exuder are activated.
  //
  // Optimization processes are designed to improve mesh quality. However, beware that such an
  // improvement is obtained by perturbing mesh vertices and modifying the mesh connectivity which
  // has an impact on the strict compliance to the refinement criteria. Though a strict compliance to
  // mesh criteria is granted at the end of the Delaunay refinement, this may no longer be true after
  // some optimization processes. Also beware that the default behavior does involve some optimization processes.
  // Any suborder of (odt -> lloyd -> perturb -> exude) is fine.
 
  // The functions lloyd_optimize_mesh_3() and odt_optimize_mesh_3()  are mesh optimization processes
  // based on the minimization of a global energy function. The minimized global energy may be
  // interpreted as the L2-norm of the error achieved when the function x^2 is interpolated on the mesh
  // domain using a piecewise linear function which is linear in each cell of the Voronoi diagram of the
  // mesh vertices.
  //
  // The optimizer lloyd_optimize_mesh_3() works in iterative steps. At each iteration, mesh vertices
  // are moved into positions that bring to zero the energy gradient and the Delaunay triangulation is
  // updated.
  // !Vertices on the mesh boundaries are handled in a special way so as to preserve an accurate
  //  representationof the domain boundaries!
  //
  // There are four optional paramater to control the process:
  // > time_limit             ... in seconds, a CPU time limit after which the optimization process is
  //                            stopped valid range:  0 <= sliver_bound <= 180
  // > max_iteration_number    ... sets a limit on the number of performed iterations
  // > convergence             ... the optimization process is stopped, when at the last iteration, the
  //                             displacement of any vertex is less than a given percentage of the
  //                            length of the shortest edge incident to that vertex. The parameter
  //                            gives the threshold ratio.
  // > freeze_bound            ... is designed to reduce running time of each optimization iteration.
  //                             Any vertex that has a displacement less than a given percentage of
  //                            the length (the of its shortest incident edge, is frozen (i.e. not
  //                            not relocated). The parameter fives the threshold ratio.
  // > do_freeze             ... If it is set to true (default value), frozen vertices will not move
  //                            anymore in next iterations. Otherwise, at each iteration, any vertex
  //                            that moves, unfreezes all its incident vertices 
  // Defaults are 'time_limit=0', 'max_iteration_number=0', meaning no early stopping, and
  // 'convergence=0.02', 'freeze_bound=0.01', 'do_freeze=true'.
  if( use_odt )
  {
      std::cout << "Performing ODT optimization... " << std::flush; 
      t_start = tic();
      CGAL::odt_optimize_mesh_3(c3t3, domain, time_limit=0, max_iteration_number=0, convergence=0.01, freeze_bound=0.005, do_freeze=true);
      t_end = tic();
      printDuration( t_start, t_end );
  }

  if( use_lloyd )
  {
      std::cout << "Performing Lloyd optimization... " << std::flush; 
      t_start = tic();
      CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=0, max_iteration_number=0, convergence=0.01, freeze_bound=0.005, do_freeze=true);
      t_end = tic();
      printDuration( t_start, t_end );
  }

  // The perturber tries to improve the dihedral angles of the worst cells in the mesh degree by degree:
  // the step number 'n' is considered as successful if after this step the worst tetrahedron of the mesh
  // has a minimal dihedral angle larger than 'n' degrees.
  // Two parameters are available to terminate the pertuber early:
  // > time_limit     ... in seconds, a CPU time limit after which the optimization process is stopped
  //                    valid range:  0 <= sliver_bound <= 180
  // > sliver_bound    ... in degree, a targeted lower bound on dihedral angles of mesh cells
  //                     (after which the worst tetrahedron in the mesh has a smallest angle larger
  //                    than 'sliver_bound degrees')
  //                    valid range: time_limit >= 0
  //
  // Default values are 'time_limit=0', 'sliver_bound=0' meaning the perturber will not stop early, but
  // run until an iteration step fails.
  // Note: in that case parameters do not need to be provided. For completeness I have provided them.
  if( perturb )
  {
      std::cout << "Pertubating mesh... " << std::flush; 
      t_start = tic();
      CGAL::perturb_mesh_3( c3t3, domain, time_limit=0, sliver_bound=0);
      t_end = tic();
      printDuration( t_start, t_end );
  }

  // The sliver exudation process consists in turning the Delaunay triangulation into a weighted
  // Delaunay triangulation and optimizing the weights of vertices in such a way that slivers disappear.
  // Two parameters are available to terminate the exuder early:
  // > time_limit     ... in seconds, a CPU time limit after which the optimization process is stopped
  //                    valid range:  0 <= sliver_bound <= 180
  // > sliver_bound ... in degree, a targeted lower bound on dihedral angles of mesh cells. The exudation
  //                    process considers in turn all the mesh ells that have a smallest dihedral angle
  //                    than 'sliver_bound' and tires to make them disappear by weighting their vertices.
  //                    The optimization process stops when every cell in the mesh achieves this quality.
  //                    
  // Default values are 'time_limit=0', 'sliver_bound=0" meaning the exuder will not stop early, but runs
  // as long as it can improve the smallest dihedral angle of the set of cells incident to some vertices.
  // Note: in that case parameters do not need to be provided. For completeness I have provided them.
  if( exude )
  {
      std::cout << "Exuding slivers... " << std::flush; 
      t_start = tic();
      CGAL::exude_mesh_3(c3t3, sliver_bound=0, time_limit=0);
      t_end = tic();
      printDuration( t_start, t_end );
  }

  // Output - Medit | MESH (temporary output)
  std::ofstream meditFile("temp.mesh");
  c3t3.output_to_medit( meditFile );
  meditFile.close();

  // Output - gmsh | MSH
  GmshInitialize( );
  GModel *m = new GModel();
  m->readMESH( "temp.mesh"  );

  // Mark all volumes as physical volumes.
  // Then they will be recognized as CellGroups in OpenFOAM. 
  int ctr = 0, electrode_ctr = 0;
  for(auto region = m->firstRegion(); region != m->lastRegion(); region++)
  {
    std::string name("Physical Volume ");
    name.append( std::to_string(ctr) );
    (*region)->addPhysicalEntity( m->setPhysicalName(name,3) );
      ctr++;
  }

  // OpenFOAM interprets Physical Surfaces as patches
  // However patches are only valid to occur on the boundary of the
  // domain, not inside. Therefore only for the electrodes a
  // definition of Physical Surfaces is valid.
  ctr = ctr - electrodepaths.size();
  //std::cout << "#physicals = " << m->numPhysicalNames() << ", #polys = " << polydomains.size() << ", counter = " << ctr << std::endl;
  for(auto face = m->firstFace(); face != m->lastFace(); face++)
  {
    
    if( ctr <= 0 )
    {
        std::string name("Electrode ");
        name.append( std::to_string( electrode_ctr ) );
        (*face)->addPhysicalEntity( m->setPhysicalName( name,2 ) );
          electrode_ctr++;
    }

    // we skip the first surfaces, because we know that the
    // electrodes are always the last surface (since we assigned
    // a subdomain index > 100 to them)
    ctr--;
  }
  m->writeMSH( outfile );
 
  // cleanup 
  delete m;
  GmshFinalize();
  std::remove("temp.msh");
  std::for_each( polydomains.begin(), polydomains.end(), [](Polyhedron_domain *pd){delete pd;});

  return EXIT_SUCCESS;
}

void printHelpText()
{
    std::cout << "This tool computes a tetrahedral 3D volume mesh based on a segmented input image.\n"
              << "Usage: ./volmesh [parameters], with the following parameters\n"
              << "\t --imagefile ... path to input segmentation image file, a labelled voxel image\n"
              << "\t --electrodefiles ... path to input electrode surface files, an OFF surface description\n"
              << "\t --outfile ... name of output file/path\n"
              << "\t --f_angbound ... minimum (inner) angle of a surface triangle\n"
              << "\t --f_sizebound ... maximum size of a surface triangle\n"
              << "\t --f_distbound ... maximum distance of triangle center from object surface\n"
              << "\t --c_reratio ... maximum ratio of a cell's radius and edge length\n"
              << "\t --c_sizebound ... maximum cell size \n"
              << "\t --lloyd ... optimize mesh using the lloyd algorithm\n"
              << "\t --odt ... optimize mesh using the odt algorithm\n"
              << "\t --perturb ... eliminate bad triangles by perturbating the mesh \n"
              << "\t --exude ... exude slivers" << std::endl;
}

bool createPolyhedronFromOFF( const char *_path, Polyhedron &_poly )
{
    // read the OFF file containing the surface description of the electrode
      std::ifstream off_input(_path, std::ifstream::in);
      if( !off_input.is_open() )
      {
        std::cout << "Cannot access surface file!" << std::endl;
          return false;
    }

    // assemble Polyhedron based on the surface description
    off_input >> _poly;
    if( off_input.fail() )
    {
        std::cerr << "Error: Cannot assemble a valid polyhedron from input file " << _path << std::endl;
        return false;
    }
    off_input.close();

    if( _poly.empty() )
    {
        std::cout << "Polyhedron is empty!" << std::endl;
        return false;
    }

    return true;
}

/* Scribbles for later use/reference of GMSH api functions:

#1 access all entities (volumes, surfaces ...)
  std::vector<GEntity*> entities;
  m->getEntities(entities); 

#2 save C3T3 as vtkUnstructuredGrid 
#include <CGAL/IO/Polyhedron_iostream.h>

  vtkUnstructuredGrid* out = CGAL::output_c3t3_to_vtk_unstructured_grid(c3t3);
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(outfile);
  writer->SetInputData(out);
  writer->Write();

#3 save surfaces as STL files ( note: not properly working)
//#include "c3t3ToSTL.hpp"

  std::ofstream offFile1("gm.off");
  std::ofstream offFile2("wm.off");
  std::ofstream offFile3("ventricle.off");
  c3t3.output_boundary_to_off( offFile1, 1 );
  c3t3.output_boundary_to_off( offFile2, 2 );
  c3t3.output_boundary_to_off( offFile3, 3 );

  SubdomainList<C3t3::Subdomain_index> subdomains;
  subdomains.push_back( SubdomainRecord<C3t3::Subdomain_index>( 1, std::string( "gm" ) )); 
  subdomains.push_back( SubdomainRecord<C3t3::Subdomain_index>( 2, std::string( "wm" ) )); 
  subdomains.push_back( SubdomainRecord<C3t3::Subdomain_index>( 3, std::string( "ventricles" ) )); 

  output_boundary_of_c3t3_to_stl<C3t3>( c3t3 , subdomains, offFile1 );

#4 save surfaces as OFF files (working)
  std::ofstream out(std::string(outfile) + std::string(".off"));
  CGAL::output_surface_facets_to_off (out, meshHolder);

*/
