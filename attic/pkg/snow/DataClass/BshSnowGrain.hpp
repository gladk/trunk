#pragma once

#include<yade/core/GeometricalModel.hpp>
#include<Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>
#include<vector>
#include<map>
#include<boost/serialization/vector.hpp>
#include<boost/serialization/shared_ptr.hpp>
#include<boost/tuple/tuple.hpp>
#include<boost/thread.hpp>


typedef std::vector< std::vector<std::vector<unsigned char> > > T_DATA;

// delete this class after we migrate to boost::serialization
class Grrrr : public Serializable
{ // a workaround to stupid bug in yade::serialization
	public:
		std::vector<Vector3r> grr;
		Grrrr(){};
		Grrrr(std::vector<Vector3r>& g):grr(g){};
	protected :
		REGISTER_ATTRIBUTES(Serializable, (grr));
	REGISTER_CLASS_NAME(Grrrr);
	REGISTER_BASE_CLASS_NAME(Serializable);
};

REGISTER_SERIALIZABLE(Grrrr);

class BshSnowGrain : public GeometricalModel
{
	public:
		boost::mutex		m_mutex;

		Vector3r center,c_axis;
		Vector3r start,end;
		Vector3r color;
		int selection;
		std::vector<std::vector<Vector3r> > slices;
		Real layer_distance;

		int m_how_many_faces;
		Real m_arbitrary_safety_coefficient,m_parallelepiped_depth;
		std::vector<boost::tuple<Vector3r,Vector3r,Vector3r,Vector3r> > m_faces; // A,B,C,normal
		std::map<Vector3r, std::set<Vector3r> > m_edges;
		std::vector<Real> m_depths; // depth for each face (allows faster checking of collision).
		// depths are negative numbers! positive number would be an altitude and means that point is _above_ the face
		size_t m_lookup_resolution;
		Vector3r m_min,m_max,m_dist;
		std::vector<std::vector<std::vector<std::set<int> > > > m_quick_lookup; // quick lookup table: x,y,z,face ids, when face ids is empty - then it's outside. /////////NOT: When face ids contain -1 - then it is inside. Otherwise it is face IDs .

		std::vector<Grrrr> gr_gr;
	public: 
		BshSnowGrain():GeometricalModel(){createIndex(); m_how_many_faces=-1;};
		BshSnowGrain(const T_DATA& dat,Vector3r c_axis,int SELECTION,Vector3r col,Real one_voxel_in_meters_is,Real layer_distance_voxels,Real angle_increment);
		Vector3r search(const T_DATA& dat,Vector3r c,Vector3r dir);
		Vector3r search_plane(const T_DATA& dat,Vector3r c,Vector3r dir);

		bool is_point_inside_polyhedron(Vector3r point);
		bool is_point_inside_polyhedron_without_quick_lookup(Vector3r P);
		bool check_if_point_is_inside_single_face(Vector3r P, size_t i);
		void check_edge_with_quick_lookup_table(Vector3r A,Vector3r B,size_t triangle_id);
		void has_deformed(){m_faces.clear();m_edges.clear();m_how_many_faces=-1;};
		int how_many_faces();
		bool face_is_valid(Vector3r&,Vector3r&,Vector3r&);
		Real depth(int i){return m_depths[i];};
		void push_face(Vector3r,Vector3r,Vector3r);
		void add_edge(Vector3r a,Vector3r b);
		const std::vector<boost::tuple<Vector3r,Vector3r,Vector3r,Vector3r> >& get_faces_const_ref(){how_many_faces(); return m_faces;};
		const std::map<Vector3r, std::set<Vector3r> >& get_edges_const_ref(){how_many_faces(); return m_edges;};
		std::vector<boost::tuple<Vector3r,Vector3r,Vector3r,Vector3r> > get_faces_copy(){how_many_faces(); boost::mutex::scoped_lock scoped_lock(m_mutex); return m_faces;};
		std::map<Vector3r, std::set<Vector3r> > get_edges_copy(){how_many_faces(); boost::mutex::scoped_lock scoped_lock(m_mutex); return m_edges;};
	
	private:
		Real calc_depth(size_t);
		bool is_point_orthogonally_projected_on_triangle(Vector3r& a,Vector3r& b,Vector3r c,Vector3r& N,Vector3r& P,Real point_plane_distance = 0.0);

		// friend class boost::serialization::access;
		template<class Archive>
		void serialize_REMOVED(Archive & ar, unsigned int version)
		{
			//boost::mutex::scoped_lock scoped_lock(m_voxel_envelope_mutex);
			ar & layer_distance;
			ar & center & c_axis;
			ar & start & end;
			ar & color;
			ar & selection;
			ar & slices;
		}
	
	protected :
		REGISTER_ATTRIBUTES(GeometricalModel,
			(center)
			(c_axis)
			(start)
			(end)
			(color)
			(selection)
			(layer_distance)
			(gr_gr) // slices
		);
		void preProcessAttributes(bool loading)
		{
			if(!loading)
			{
				gr_gr.clear();
				BOOST_FOREACH(std::vector<Vector3r>& g,slices)
					gr_gr.push_back(Grrrr(g));
			}
		};
		void postProcessAttributes(bool loading)
		{
			if(loading)
			{
				m_how_many_faces = -1;
				slices.clear();
				BOOST_FOREACH(Grrrr& g,gr_gr)
					slices.push_back(g.grr);
			}
		};
	REGISTER_CLASS_NAME(BshSnowGrain);
	REGISTER_BASE_CLASS_NAME(GeometricalModel);
	REGISTER_CLASS_INDEX(BshSnowGrain,GeometricalModel);
};

REGISTER_SERIALIZABLE(BshSnowGrain);

