#include"BshSnowGrain.hpp"
#include<Wm3Quaternion.h>

YADE_REQUIRE_FEATURE(geometricalmodel);


// a voxel is 20.4 microns (2.04 × 10-5 meters)
// the sample was 10.4mm hight


BshSnowGrain::BshSnowGrain(const T_DATA& dat,Vector3r c_ax,int SELECTION,Vector3r col, Real one_voxel_in_meters_is,Real layer_distance_voxels,Real angle_increment) : GeometricalModel()
{
	createIndex();

	layer_distance = layer_distance_voxels;

	if(SELECTION!=0)
	{
		color=col;

		c_axis=c_ax;
		c_axis.Normalize();
		selection=SELECTION;
		std::cout << "creating grain " << SELECTION << " c_axis " << c_axis << "\n";

		Vector3r sum(0,0,0);
		float count=0;

		int FRAMES=dat.size();
		int sx=dat[0].size();
		int sy=dat[0][0].size();
		for(int i=0;i<FRAMES;++i)
		{
			for(int x=0;x<sx;++x)
				for(int y=0;y<sy;++y)
					if( dat[i][x][y]==SELECTION )
					{
						sum+=Vector3r((float)x,(float)y,(float)i);
						count+=1.0;
					}
		}
		center=sum/count;

		std::cout << "finding start/end\n";

		start=search_plane(dat,center,c_axis);
		end  =search_plane(dat,center,-1.0*c_axis);
		
		std::cout << "(start-end).length() " << ((start-end).Length()) << " layer_distance (voxels): " << layer_distance << "\n";

		Quaternionr q;
		q.Align(Vector3r(0,0,1),c_axis);
//////////////////////////////////////////////////////////
int II=0;
float prevL=(start-end).SquaredLength();
float tmpL;
Vector3r moving_center(0,0,0),moving_sum(0,0,0);
count=0;
//layer_distance=3.0;
Vector3r S=start;
std::cerr << "\n--------- " << (int)((bool)((tmpL=(S-end).SquaredLength())<=prevL)) << "\n";
std::cerr << tmpL << " " << prevL << "\n";
for(; (tmpL=(S-end).SquaredLength())<=prevL ; S-=c_axis*layer_distance, ++II)
{
	slices.resize(II+1);
	float vectorY1=5.0;
	float vectorX1=0.0;
	for(float angle=0.0f ; angle <= (2.0f*3.14159) ; angle+=angle_increment)
	{		
		float vectorX=(5.0*sin(angle));
		float vectorY=(5.0*cos(angle));		
		Vector3r tmp(vectorX1,vectorY1,0);
		tmp=search(dat,S+moving_center,q*tmp);
		slices[II].push_back(tmp-center);
		vectorY1=vectorY;
		vectorX1=vectorX;

		moving_sum+=tmp;
		count+=1.0;
	}
	prevL=tmpL;

	moving_center=moving_sum/count-S;
	moving_sum=Vector3r(0,0,0);
	count=0;

	std::cout << "slices " << II << " S-end " << (S-end).SquaredLength() << "\n";
}
//////////////////////////////////////////////////////////
// we're done. Now scale everything.
BOOST_FOREACH(std::vector<Vector3r>& g,slices)
	BOOST_FOREACH(Vector3r& v,g)
		v *= one_voxel_in_meters_is;

c_axis.Normalize();
start = (start-center)*one_voxel_in_meters_is;
end = (end-center)*one_voxel_in_meters_is;
center *= one_voxel_in_meters_is;


		
		std::cout << "grain created, center " << center << " start " << start << " end " << end << "\n";
	} 
	else
		std::cout << "nothing selected!\n";
};
	
Vector3r BshSnowGrain::search(const T_DATA& dat,Vector3r c,Vector3r dir)
{
	dir.Normalize();
	dir=dir*0.3;

	Vector3r the_farthest(c);

	int FRAMES=dat.size();
	int sx=dat[0].size();
	int sy=dat[0][0].size();
	for( Vector3r current(c); 
		   ((int)(current[0])<sx    ) && ((int)(current[0])>=0) 
		&& ((int)(current[1])<sy    ) && ((int)(current[1])>=0) 
		&& ((int)(current[2])<FRAMES) && ((int)(current[2])>=0) 
		; current+=dir )

			if( dat[(int)(current[2])][(int)(current[0])][(int)(current[1])]==selection )
				the_farthest=current;

	return the_farthest;
};

Vector3r BshSnowGrain::search_plane(const T_DATA& dat,Vector3r c,Vector3r dir)
{
	Vector3r result=search(dat,c,dir);
	
	dir.Normalize();
	
	Quaternionr q;
	q.Align(Vector3r(0,0,1),dir);

	dir=dir*0.3;
	bool anything_found=true;
		
	for( ; anything_found==true ; result+=dir)
	{
		anything_found=false;

		float vectorY1=5.0;
		float vectorX1=0.0;
		for(float angle=0.0f ; (angle <= (2.0f*3.14159)) && (anything_found==false) ; angle+=0.2f)
		{		
			float vectorX=(5.0*sin(angle));
			float vectorY=(5.0*cos(angle));		
			Vector3r tmp(vectorX1,vectorY1,0);
			tmp=search(dat,result,q*tmp);

			if((tmp-result).SquaredLength()>1.0)
				anything_found=true;
			else
				anything_found=false;

			vectorY1=vectorY;
			vectorX1=vectorX;			
		}

	}

	return result;
};

bool BshSnowGrain::is_point_orthogonally_projected_on_triangle(Vector3r& a,Vector3r& b,Vector3r c,Vector3r& N,Vector3r& P,Real point_plane_distance)
{
	// first check point - plane distance
	if(point_plane_distance == 0.0)
		point_plane_distance = N.Dot(P - c);
	if( point_plane_distance <= 0 )            // point has to be inside - on negative side of the plane 
	{
		// now calculate projection of point in the plane
		Vector3r d(P - point_plane_distance*N);
		// now check if the point (when projected on a plane) is within triangle a,b,c
		// it could be faster with methods from http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm
		// but I don't understand them, so I prefer to use the method which I derived myself
		Vector3r c1((a - b).Cross(d - a)); // since points a,b,c are all clockwise, 
		Vector3r c2((c - a).Cross(d - c)); // then if I put a point 'd' inside a triangle it will be clockwise
		Vector3r c3((b - c).Cross(d - b)); // with each pair of other points, but will be counterclockwise if it's outside
		if(c1.Dot(N) >= 0 && c2.Dot(N) >= 0 && c3.Dot(N) >= 0) // therefore if any of them is counterclockwise, the dot product will be negative
			return true;
	}
	return false;
};

bool BshSnowGrain::check_if_point_is_inside_single_face(Vector3r P, size_t i)
{
	Vector3r a(get<0>(m_faces[i]));
	Vector3r b(get<1>(m_faces[i]));
	Vector3r c(get<2>(m_faces[i]));
	Vector3r N(get<3>(m_faces[i]));
	Real depth = m_depths[i];
	Real point_plane_distance = N.Dot(P - c);
	if(   point_plane_distance < 0             // point has to be inside - on negative side of the plane
	   && point_plane_distance > depth ) // and has to be within the depth of this face 
	{
		if(is_point_orthogonally_projected_on_triangle(a,b,c,N,P,point_plane_distance))
		{
			if( point_plane_distance > depth*m_parallelepiped_depth ) // (1-parallelepiped)
			{ // that's close enough. We can speed up the computation by returning true at this point

// therefore in fact, we are checking if point is inside a volume of parallelepiped (1-) the height of 1/2 depth
// within the polyhedron PLUS (that's the code below - the (2-)) a polyhedron the height given by depth
// like this:
//              ____________     this weird shape is WHOLE polyhedron
//             /         Z. \                                                         ------
//            /          ' ` \                                                             |
//           /          '2-tetrahedron with the height equal to depth of THIS trangle      |
//          /          '     ` \                                                           | 
//         /          '       ` \                                                          |d
//        /          '         ` \______________ XX face (see description)                 |e - search below for
//       /          '           `              /                                           |p 'arbitrary safety coefficient'
//      /      ....'.............`......      /                                            |t
//      \     '   '           1-parallelepiped with height equal to 1/2 of depth           |h ('Z' is orthogonal to triangle with 
//       \    '  '                 `   '    /                                              |   'depth' distance from the
//        \   ' '                   `  '   /                                               |   triangle's center point )
//         \  ''                     ` '  /                                                |
//          \ '                       `' /                                                 |
//           \__________________________/                                             ------
//               THIS triangular face                
//
// depth of a triangle is the distance from it to the opposite side of polyhedron, kind of
// diameter in a sphere. The depth is calculated in four places: at three nodes of a triangle,
// and in its center point. The shallowest one is chosen as the actual value. This is
// not describing volume in exactly perfect math, but good enough, and rather fast.
// Sometimes for example the 'depth' of some triangle can extend the volume of the real polyhedron a bit
// on the other side. For example if that 'XX face' was closer to 'THIS triangular face' it would 
// still be considered as inside of polyhedron, because the '1/2 depth parallelepiped' might extend beyond it.
// this could happen if none of four (points a,b,c and center point) reached 'XX face' when calculating depth,
// because it was off a bit, and depth calculation ended up on another face.

// From this approach you can see that it's optimized for polyhedrons with large number of faces. For instance it will
// be terribly wrong if polyhedron is just a single four-noded tetrahedron. But it will work well if polyhedron
// is made from many triangles.
				return true;
			}
			else
			{ // (2-terhahdron)
				// so a point orthogonally projected on triangle is crossing it
				// now check if this point is inside a tetrahedron made by this triangle
				// and a point 'Z' at 'depth', see picture
				Vector3r Z((a+b+c)/3.0 + N*depth);
				Vector3r N1((Z - a).Cross(a - b)); // normals of each face
				Vector3r N2((Z - c).Cross(c - a)); // of tetrahedron a,b,c,Z
				Vector3r N3((Z - b).Cross(b - c));
				if(
					is_point_orthogonally_projected_on_triangle(b,a,Z,N1,P) &&
					is_point_orthogonally_projected_on_triangle(a,c,Z,N2,P) &&
					is_point_orthogonally_projected_on_triangle(c,b,Z,N3,P)
					)
					return true;
			}
		}
	}
	return false;
}

bool BshSnowGrain::is_point_inside_polyhedron_without_quick_lookup(Vector3r P)
{
	how_many_faces();
	// loop on all faces
	size_t S(m_faces.size());
	for(size_t i = 0; i < S ; ++i)
		if(check_if_point_is_inside_single_face(P,i))
			return true;
	return false;
};

bool BshSnowGrain::is_point_inside_polyhedron(Vector3r P)
{
//	return is_point_inside_polyhedron_without_quick_lookup(P);

	how_many_faces();
	Vector3r Q(P-m_min);
	int i = std::floor(Q[0]/m_dist[0]);
	int j = std::floor(Q[1]/m_dist[1]);
	int k = std::floor(Q[2]/m_dist[2]);
	if(i<0 || j<0 || k<0 || i>= (int)m_lookup_resolution || j>= (int)m_lookup_resolution || k>= (int)m_lookup_resolution)
		return false;
	std::set<int>& ids(m_quick_lookup[i][j][k]);
	if(ids.empty())
		return false;
//	if(ids.size()==1 && *(ids.begin())==-1)
//		return true;
	BOOST_FOREACH(int i,ids)
		if(/*i != -1 && */check_if_point_is_inside_single_face(P,i))
			return true;
	return false;
};

bool BshSnowGrain::face_is_valid(Vector3r& a,Vector3r& b,Vector3r& c)
{
	if(a != b && b != c && c != a)
		return true;
	return false;
};

void BshSnowGrain::add_edge(Vector3r a,Vector3r b)
{
	if(m_edges.find(a) == m_edges.end())
		m_edges[a]=std::set<Vector3r>();
	if(m_edges.find(b) == m_edges.end())
		m_edges[b]=std::set<Vector3r>();
	m_edges[a].insert(b);
	m_edges[b].insert(a);
/*	
	typedef std::pair<Vector3r,Vector3r> t_edge;
//	std::map<Vector3r, std::set<Vector3r> > m_edges;
	t_edge ab(std::make_pair(a,b));
	t_edge ba(std::make_pair(b,a));

	int a_first(0),b_first(0);
	BOOST_FOREACH(t_edge& e,m_edges)
	{
		if(e == ab || e == ba) return;
		if(e.first == a) ++a_first;
		if(e.first == b) ++b_first;
	}
	if((a_first == 0 && b_first == 0) || (a_first == 0 && b_first == 1))
	{
		m_edges.push_back(ab);
		return;
	}
	if(a_first == 1 && b_first == 0)
	{
		m_edges.push_back(ba);
		return;
	}
	if(a_first == 1 && b_first == 1)
	{
		std::cerr << "BshSnowGrain::add_edge error, both nodes already first.\n";
		return;
	}
	std::cerr << "BshSnowGrain::add_edge error, both nodes present multiple times!!.\n";
*/
}

void BshSnowGrain::push_face(Vector3r a,Vector3r b,Vector3r c)
{
	if(face_is_valid(a,b,c))
	{
		Vector3r n((a - b).Cross(c - a));
		if(n.SquaredLength() != 0)
		{
			n /= n.Length();
			m_faces.push_back(boost::make_tuple(a,b,c,n));
			add_edge(a,b);
			add_edge(b,c);
			add_edge(c,a);
		} else
		{
			std::cerr << "Face has no normal!\n";
			n=Vector3r(1,0,0);
		}
	}
}
		
int BshSnowGrain::how_many_faces()
{
	//boost::mutex::scoped_try_lock scoped_try_lock(m_mutex,boost::try_to_lock_t);
	//if(! scoped_try_lock.owns_lock())
	//{
	//	std::cerr << "..lock failed\n";
	//}

	if(m_how_many_faces != -1)
		return m_how_many_faces;

// sorry, I never got around to find time and to check out how LOG_WARN works .... std::cerr is fine for me ;)
//	std::cerr << "recalculating the depths of polyhedron triangular faces - for faster collision detection\n";


// CREATE TRIANGULAR FACES. usually a polyhedron has triangular faces, but here it's a snow grain.
// it has layers, not faces, I have to make faces from layers

	m_faces.clear();
	m_edges.clear();
	//calculate amount of faces..

	// connected to START - the middle point in first layer
	int S = slices[0].size();
	Vector3r START(slices[0][0]);
	for(int j = 1 ; j < S ; ++j)
		START += slices[0][j];
	START /= (float)(S);
	for(int j = 0 ; j < S ; ++j)
		push_face( slices[0][j] , slices[0][(j+1 < S) ? (j+1):0] , START );

	// all triangles between layers
	int L = slices.size();
	for(int i = 1 ; i < L ; ++i)
		for(int j = 0 ; j < S/*slices[i].size()*/ ; ++j)
		{
			push_face( slices[i][j] , slices[i-1][j] , slices[i-1][(j-1 > 0) ? (j-1):(S-1)] );
			push_face( slices[i][j] , slices[i][(j+1 < S) ? (j+1):0] , slices[i-1][j] );
		}
	
	// connected to END - the middle point in last layer
	Vector3r END(slices[L-1][0]);
	for(int j = 1 ; j < S ; ++j)
		END += slices[L-1][j];
	END /= (float)(S);
	for(int j = 0 ; j < S ; ++j)
		push_face( slices[L-1][j] , END , slices[L-1][(j+1 < S) ? (j+1):0] );
	
	m_how_many_faces = m_faces.size();


// NOW I HAVE FACES. That's the code for usual polyhedrons, that already have faces:
// calculating the depth for each face.

	// now calculate the depth for each face
	m_arbitrary_safety_coefficient = 0.7;
	m_parallelepiped_depth = 0.5;
	m_depths.resize(m_faces.size(),0);
	// loop on all faces
	size_t SS(m_faces.size());
	for(size_t i = 0; i < SS ; ++i)
		m_depths[i] = calc_depth(i)*m_arbitrary_safety_coefficient; // 0.7 is an arbitrary safety coefficient

// NOW let's calculate quick spatial lookup table

	m_lookup_resolution=15;
	// resize the table.
	m_quick_lookup.resize(m_lookup_resolution);
	BOOST_FOREACH(std::vector<std::vector<std::set<int> > >& a,m_quick_lookup)
	{
		a.resize(m_lookup_resolution);
		BOOST_FOREACH(std::vector<std::set<int> >& b,a) b.resize(m_lookup_resolution);
	}
	// calculate min/max of Aabb of polyhedron
	m_min = get<0>(m_faces[0]);
	m_max = m_min;
	for(size_t i = 0; i < SS ; ++i)
	{
		Vector3r A(get<0>(m_faces[i]));
		Vector3r B(get<1>(m_faces[i]));
		Vector3r C(get<2>(m_faces[i]));
	 	m_max = componentMaxVector(m_max,A);m_min = componentMinVector(m_min,A);
	 	m_max = componentMaxVector(m_max,B);m_min = componentMinVector(m_min,B);
	 	m_max = componentMaxVector(m_max,C);m_min = componentMinVector(m_min,C);
	}
	m_dist=(m_max-m_min)/((Real)(m_lookup_resolution));
//	// first a rough check - mark whatever is inside.
//	for(size_t i=0 ; i< m_lookup_resolution-1 ; ++i)
//	for(size_t j=0 ; j< m_lookup_resolution-1 ; ++j)
//	for(size_t k=0 ; k< m_lookup_resolution-1 ; ++k)
//	{
//		Vector3r pos_min(i  ,j  ,k  );pos_min=diagMult(pos_min,m_dist)+m_min;
//		Vector3r pos_max(i+1,j+1,k+1);pos_max=diagMult(pos_max,m_dist)+m_min;
//		Vector3r pos_1(pos_min[0],pos_min[1],pos_max[2]);
//		Vector3r pos_2(pos_min[0],pos_max[1],pos_max[2]);
//		Vector3r pos_3(pos_min[0],pos_max[1],pos_min[2]);
//		Vector3r pos_4(pos_max[0],pos_min[1],pos_max[2]);
//		Vector3r pos_5(pos_max[0],pos_min[1],pos_min[2]);
//		Vector3r pos_6(pos_max[0],pos_max[1],pos_min[2]);
//		if( // check all 8 corners
//			   is_point_inside_polyhedron_without_quick_lookup(pos_min)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_max)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_1)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_2)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_3)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_4)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_5)
//			&& is_point_inside_polyhedron_without_quick_lookup(pos_6)
//		)
//			m_quick_lookup[i][j][k].insert(-1);
//	}
	// now check all edges of the "triangle"'s volume
//	Real D = m_dist.Length();
	for(size_t i = 0; i < SS ; ++i)
	{
		Vector3r a(get<0>(m_faces[i]));
		Vector3r b(get<1>(m_faces[i]));
		Vector3r c(get<2>(m_faces[i]));
		Vector3r n(get<3>(m_faces[i]));
		Real depth = m_depths[i];
		// face of triangle
		check_edge_with_quick_lookup_table(a,b,i);
		check_edge_with_quick_lookup_table(b,c,i);
		check_edge_with_quick_lookup_table(c,a,i);
		
		Vector3r Z((a+b+c)/3.0 + n*depth);
		// its internal tetrahedron
		check_edge_with_quick_lookup_table(a,Z,i);
		check_edge_with_quick_lookup_table(b,Z,i);
		check_edge_with_quick_lookup_table(c,Z,i);
		// its internal parallelepiped
		Real depth2 = depth*m_parallelepiped_depth;
		Vector3r N(n*depth2);
		Vector3r A(a+N);
		Vector3r B(b+N);
		Vector3r C(c+N);
		// its top base face
		check_edge_with_quick_lookup_table(A,B,i);
		check_edge_with_quick_lookup_table(B,C,i);
		check_edge_with_quick_lookup_table(C,A,i);
		
		// all connections between its nodes
		check_edge_with_quick_lookup_table(A,a,i);
		//check_edge_with_quick_lookup_table(A,b,i);
		//check_edge_with_quick_lookup_table(A,c,i);
		
		//check_edge_with_quick_lookup_table(B,a,i);
		check_edge_with_quick_lookup_table(B,b,i);
		//check_edge_with_quick_lookup_table(B,c,i);
		
		//check_edge_with_quick_lookup_table(C,a,i);
		//check_edge_with_quick_lookup_table(C,b,i);
		check_edge_with_quick_lookup_table(C,c,i);
	}

	return m_how_many_faces;
};

void BshSnowGrain::check_edge_with_quick_lookup_table(Vector3r A,Vector3r B,size_t triangle_id)
{
	Vector3r Z(A);
	int sections = m_lookup_resolution;
	Vector3r DD((B-A)/((Real)(sections)));
	for(int zzz=0 ; zzz<=sections ; ++zzz,Z+=DD)
	{
		Vector3r Q(Z-m_min);
		int i = std::floor(Q[0]/m_dist[0]);
		int j = std::floor(Q[1]/m_dist[1]);
		int k = std::floor(Q[2]/m_dist[2]);
		if(i<0) i=0;
		if(j<0) j=0;
		if(k<0) k=0;
		if(i>=(int)m_lookup_resolution) i=m_lookup_resolution-1;
		if(j>=(int)m_lookup_resolution) j=m_lookup_resolution-1;
		if(k>=(int)m_lookup_resolution) k=m_lookup_resolution-1;
		m_quick_lookup[i][j][k].insert(triangle_id);
	}
}

Real BshSnowGrain::calc_depth(size_t I)
{
	Vector3r A(get<0>(m_faces[I]));
	Vector3r B(get<1>(m_faces[I]));
	Vector3r C(get<2>(m_faces[I]));
	Vector3r N(get<3>(m_faces[I]));
	Vector3r P((A+B+C)/3.0);
	// (1) ray N is cast from point P, where P is on some triangle.
	// return the distance from P to next closest triangle

	Real depth = 0;
	const std::vector<boost::tuple<Vector3r,Vector3r,Vector3r,Vector3r> >& f(get_faces_const_ref());
	// loop on all faces
	size_t S(f.size());
	for(size_t i = 0; i < S ; ++i)
	{
		if(I != i) // don't check with itself
		{
			Vector3r n(get<3>(f[i]));
			Real parallel = n.Dot(N); // 'N' parallel to 'n' gives 0 dot product
			if( parallel < 0) // (2) must face in opposite directions
			{
				Vector3r a(get<0>(f[i]));
				Vector3r b(get<1>(f[i]));
				Vector3r c(get<2>(f[i]));
				for(int Z = 0 ; Z < 3 ; ++Z) // (ad. 1) OK, in fact it's not just from 'P' - we cast four rays.
				                             // From 'P' and all triangle nodes
				{
					Vector3r PP;
					switch(Z)
					{
					//	case 0 : PP = P; break;
					//	case 1 : PP = A; break;
					//	case 2 : PP = B; break;
					//	case 3 : PP = C; break;

						case 0 : PP = (P+A)*0.5; break;
						case 1 : PP = (P+B)*0.5; break;
						case 2 : PP = (P+C)*0.5; break;
					}
					Real neg_point_plane_distance = n.Dot(c - PP);
					if( neg_point_plane_distance > 0 ) // (ad. 2) must be facing towards each other
					{
						// now calculate intersection point 'd' of ray 'N' from point 'PP' with the plane
						Real u = neg_point_plane_distance/parallel;
						Vector3r d(PP + u*N);
						// now check if the point 'd' (when projected on a plane) is within triangle a,b,c
						Vector3r c1((a - b).Cross(d - a));
						Vector3r c2((c - a).Cross(d - c));
						Vector3r c3((b - c).Cross(d - b));
						if(c1.Dot(n) > 0 && c2.Dot(n) > 0 && c3.Dot(n) > 0)
						{
							if(depth == 0)
							{ 
								depth = u;
							} else {
								depth = std::max(depth , u ); // get the shallowest one
							}
						}
					}
				}
			}
		}
	}
	return depth;
}

YADE_PLUGIN((BshSnowGrain)(Grrrr));

YADE_REQUIRE_FEATURE(PHYSPAR);

