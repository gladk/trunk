#pragma once

#include "Config.hpp"
#include "DataVoxel.hpp"

#include <vector>
#include <list>

// coordinates i,x,y - must allow the size of "all voxels cube"
//            i/x/y
typedef unsigned short							t_coord; 

//           (unique)       i/x/y   x/y/i
typedef std::set<std::pair<t_coord,t_coord> >				t_projection;

struct vc
{
		t_coord		level;
		t_voxel		color; // if different than m_own_id -> it's a contact
		vc(t_coord l,t_voxel c) : level(l), color(c) {};

	private:
		vc() {};
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & level;
			ar & color;
		}
};

typedef std::vector< std::list< vc > >					t_vis_slice;
//      i/x/y+min   x/y/i+min
typedef std::vector< t_vis_slice >					t_visibility;

enum enum_side {I_TOP,I_BOT,X_TOP,X_BOT,Y_TOP,Y_BOT};

class GrainSurface
{
	private:
		t_voxel		m_own_id;

		t_projection	m_on_i, // x,y  -->  i
				m_on_x, // y,i  -->  x
				m_on_y; // i,x  -->  y

		t_coord		m_min_i,
				m_min_x,
				m_min_y;

		t_coord		m_max_i,
				m_max_x,
				m_max_y;

		t_visibility	m_vis_i_top, // x,y  -->  i
				m_vis_x_top, // y,i  -->  x
				m_vis_y_top; // i,x  -->  y
		t_visibility	m_vis_i_bot, // x,y  -->  i
				m_vis_x_bot, // y,i  -->  x
				m_vis_y_bot; // i,x  -->  y

	private:
		GrainSurface() {};
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & m_own_id;
			ar & m_on_i;
			ar & m_on_x;
			ar & m_on_y;
			ar & m_min_i;
			ar & m_min_x;
			ar & m_min_y;
			ar & m_max_i;
			ar & m_max_x;
			ar & m_max_y;
			ar & m_vis_i_top;
			ar & m_vis_x_top;
			ar & m_vis_y_top;
			ar & m_vis_i_bot;
			ar & m_vis_x_bot;
			ar & m_vis_y_bot;
		}
	public:
		GrainSurface(t_voxel_data& dat, t_voxel grain_id);

		void		draw(SafeVectors3& colors,Config& c);
		void		draw_single_side(t_coord i,t_coord x,t_coord y,enum_side side,float d2);
		void		plot_sides(t_visibility& vis,t_coord A,t_coord B,t_coord min_A,t_coord min_B,enum_side side,Config& c,SafeVectors3& colors);

		t_voxel		id(){return m_own_id;};
};

