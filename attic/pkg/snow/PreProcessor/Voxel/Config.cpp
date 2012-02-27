YADE_REQUIRE_FEATURE(COMPILED_BY_INCLUDE_FROM_SnowVoxelsLoader.cpp)
#include "Config.hpp"

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include <iostream>
namespace po = boost::program_options;

Config::Config(int argc, char** argv) : 
	m_draw_with_names(false),
	m_display_selection(false),
	m_time_draw_limit(true),
	m_draw_using_surfaces(false),
	m_draw_backfaces(false),
	m_draw_surface_or_contact_or_all(2),
	m_auto_center(false),
	m_transparency(false),
	opts("Snow data reader options")
{
	opts.add_options()
		("help,h","display this help.")
		(",a",po::value<std::string>(&m_c_axis_file)->default_value("axec2/c_axis_sam8_dct1.1.txt"),
			"c-axis data file, each line in the file has: ID X Y Z.")
		(",c",po::value<std::string>(&m_color_file)->default_value("axec2/c_axismap_sam8_dct1.1.txt"),
			"colors data file, each line corresponds to subsequent ID, first column is ignored: ignored X Y Z.")
		(",f",po::value<std::string>(&m_file_basename)->default_value("txt2/sam8_dct1_250708_"),
			"basename of .txt data file to be loaded.")
		(",s",po::value<int>(&m_start)->default_value(1),
			"start slice number to be loaded.")
		(",e",po::value<int>(&m_end)->default_value(610),
			"end slice number to be loaded.")
		(",n",po::value<unsigned int>(&m_digits)->default_value(4),
			"number of digits in slice name eg. four is: _0000.txt.")
		(",g",po::value<int>(&m_grains)->default_value(120),
			"number of grains the sample.")
		(",d",po::value<int>(&m_display_steps)->default_value(12),
			"drawing display steps, for 1 - draw each voxel.")
		(",L",po::value<std::string>(&m_load_file)->default_value(""),
			"Whether to load some file on start.")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);

	if (vm.count("help")) {
            std::cout << opts << "\n";
            exit(0);
        }
}

void Config::light_selection(int grain)
{
	GLfloat matAmbient[4];
	if(selected(grain))
	{
		matAmbient[0] = 0.2;
		matAmbient[1] = 0.2;
		matAmbient[2] = 0.2;
		matAmbient[3] = 0.0;

	} else
	{
		matAmbient[0] = 0.0;
		matAmbient[1] = 0.0;
		matAmbient[2] = 0.0;
		matAmbient[3] = 0.0;
	}
	glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,matAmbient);
}

void Config::remove_transparency()
{
	if(m_transparency)
	{
		glDisable(GL_BLEND);
//		glEnable(GL_DEPTH_TEST);
	}
	light_selection(-1);
}

void Config::apply_transparency()
{
	if(m_transparency)
	{
	//	if(m_draw_backfaces)
	//	else
		if(m_draw_using_surfaces)
			glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		else
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
//		glDisable(GL_DEPTH_TEST);
	}
	else
	{
		glDisable(GL_BLEND);
//		glEnable(GL_DEPTH_TEST);
	}
//                                                                                               backfaces      no-backfaces
//               src                             dest                           alpha:         0.8  0.5  0.2   0.8  0.5  0.2
//glBlendFunc(GL_ONE                    ,     GL_SRC_ALPHA                         );break;//   ?+  ++   --     ?+  ++   ++     1
//glBlendFunc(GL_ONE                    ,     GL_ONE_MINUS_SRC_ALPHA               );break;//   --  ++   ++     ++  ++   -+     1
//
//glBlendFunc(GL_SRC_ALPHA              ,     GL_SRC_ALPHA                         );break;//   ?+  -?   --     ?+  -?   --     2
//glBlendFunc(GL_SRC_ALPHA              ,     GL_ONE_MINUS_SRC_ALPHA               );break;//   --  -?   ??     ++  -?   ??     2
//glBlendFunc(GL_ONE_MINUS_SRC_ALPHA    ,     GL_SRC_ALPHA                         );break;//   --  -?   --     +?  -?   ++     2
//
//glBlendFunc(GL_DST_ALPHA              ,     GL_SRC_ALPHA                         );break;//   ?+  ++   --     ?+  ++   ++     3
//glBlendFunc(GL_DST_ALPHA          *   ,     GL_ONE_MINUS_SRC_ALPHA               );break;//   --  ++   ++     ++  ++   ?+     3

// suggested by nam page             GL_SRC_ALPHA              ,     GL_ONE_MINUS_SRC_ALPHA
}


YADE_REQUIRE_FEATURE(PHYSPAR);

