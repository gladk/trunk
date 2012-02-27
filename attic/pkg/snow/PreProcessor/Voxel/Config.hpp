#pragma once

#include<yade/lib-QGLViewer/qglviewer.h>
#include<boost/program_options.hpp>
#include<boost/foreach.hpp>
#include<string>
#include<set>
#include<boost/archive/binary_oarchive.hpp>
#include<boost/archive/binary_iarchive.hpp>
#include<boost/serialization/vector.hpp>
#include<boost/serialization/set.hpp>
#include<boost/serialization/list.hpp>
#include<boost/serialization/utility.hpp>

#include<boost/iostreams/filtering_stream.hpp>
#include<boost/iostreams/filter/gzip.hpp>
#include<boost/iostreams/filter/bzip2.hpp>
#include<boost/iostreams/device/file.hpp>


class Config
{
	private:
		std::string	m_c_axis_file;
		std::string	m_color_file;
		std::string	m_file_basename;
		int		m_start;
		int		m_end;
		unsigned int	m_digits;
		int		m_grains;
		int		m_display_steps;
		std::set<int>	m_selection;
		bool		m_draw_with_names;
		bool		m_display_selection;
		bool		m_time_draw_limit;
		bool		m_draw_using_surfaces;
		qglviewer::Vec	m_view_direction;
		bool		m_draw_backfaces;
		signed char	m_draw_surface_or_contact_or_all;
		bool		m_auto_center;
		int		m_transparency;

		std::string	m_load_file;
		boost::program_options::options_description opts;

	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & m_c_axis_file;
			ar & m_color_file;
			ar & m_file_basename;
			ar & m_start;
			ar & m_end;
			ar & m_digits;
			ar & m_grains;
			ar & m_display_steps;
			ar & m_selection;
			ar & m_draw_with_names;
			ar & m_display_selection;
			ar & m_time_draw_limit;
			ar & m_draw_using_surfaces;
			ar & m_view_direction;
			ar & m_draw_backfaces;
			ar & m_draw_surface_or_contact_or_all;
			ar & m_auto_center;
			ar & m_transparency;
		}
	public:
		Config(int argc, char** argv);

		std::string	axis_file(){return m_c_axis_file;};
		std::string	color_file(){return m_color_file;};
		std::string	file_basename(){return m_file_basename;};
		int		start(){return m_start;};
		int		end(){return m_end;};
		unsigned int	digits(){return m_digits;};
		int		grains(){return m_grains;};

		int		display_steps(){return m_display_steps;};
		void		display_steps(int d){m_display_steps=d;};

		bool		selected(){return m_selection.size()!=0;};
		bool		selected(int s){return m_selection.find(s)!=m_selection.end();};
		void		select(int s){m_selection.insert(s);};
		void		unselect(int s){m_selection.erase(s);};
		void		clear_selection(){m_selection.clear();};
		std::set<int>	selection_copy(){return m_selection;};
		size_t		selection_size(){return m_selection.size();};

		bool		draw_with_names(){return m_draw_with_names;};
		void		draw_with_names(bool b){m_draw_with_names=b;};
		void		toggle_display_selection(){m_display_selection=!m_display_selection;};
		bool		display_selection(){return m_display_selection;};
		void		light_selection(int grain);

		void		time_draw_limit(bool t){m_time_draw_limit=t;};
		bool		time_draw_limit(){return m_time_draw_limit;};

		void		draw_using_surfaces(bool t){m_draw_using_surfaces=t;};
		bool		draw_using_surfaces(){return m_draw_using_surfaces;};

		void		view_direction(qglviewer::Vec c){m_view_direction=c;};
		qglviewer::Vec	view_direction(){return m_view_direction;};

		void		draw_backfaces(bool t){m_draw_backfaces=t;};
		bool		draw_backfaces(){return m_draw_backfaces;};

		void		toggle_draw_surface_or_contact_or_all(){ m_draw_surface_or_contact_or_all = ((m_draw_surface_or_contact_or_all+1)%3);};
		short int	draw_surface_or_contact_or_all(){return m_draw_surface_or_contact_or_all;};

		void		toggle_auto_center(){m_auto_center = !m_auto_center;};
		bool		auto_center(){return m_auto_center;};

		void		toggle_transparency(){m_transparency = !m_transparency;};
		int		transparency(){return m_transparency;};
		void		apply_transparency();
		void		remove_transparency();

		std::string	load_file(){return m_load_file;};
};

