YADE_REQUIRE_FEATURE(COMPILED_BY_INCLUDE_FROM_SnowVoxelsLoader.cpp)
#include <boost/bind.hpp>
#include "VoxelEnvelope.hpp"

VoxelEnvelope::VoxelEnvelope(Config& c, std::string name) : 
	m_centers_calculated(false),
	m_grain_voxel_surfaces_calculated(false),
	m_bad_grain_count(""),
	m_myname(name),
	m_data_voxel(), 
	m_axes(c.grains(),"c axes",c.axis_file(),true), 
	m_colors(c.grains(),"colors",c.color_file(),false)
{
	m_centers.clear();
}

void VoxelEnvelope::calculate_surfaces_loop(Config& c)
{
	do {
		calculate_grain_voxel_surfaces(c);
	} while(!grain_voxel_surfaces_calculated());
	std::cerr << "\nsurfaces calculated\n";
}

void VoxelEnvelope::calculate_surfaces_loop_spawn(Config& c)
{
	boost::thread loop(boost::bind(&VoxelEnvelope::calculate_surfaces_loop,boost::ref(*this),boost::ref(c)));
#if BOOST_VERSION >= 103500
	loop.detach();
#endif
}

void VoxelEnvelope::save_txt(Config& c)
{
	m_data_voxel.save_txt_result(c);
}

void VoxelEnvelope::load_experimental_data(Config& c)
{
	m_myname = "experimantal data";
	m_axes	= SafeVectors3(c.grains(),"c axes",c.axis_file(),true); 
	m_colors = SafeVectors3(c.grains(),"colors",c.color_file(),false);
	m_centers_calculated = false;
	m_grain_voxel_surfaces_calculated = false;
	m_bad_grain_count = "";
	m_centers.clear();

//	calculate_surfaces_loop_spawn(c);
	m_data_voxel.load(c);
	calculate_grain_voxel_surfaces(c);
}

void VoxelEnvelope::create_singled_particles(VoxelEnvelope& other,Config& config,std::string name)
{
	m_myname = name;
	m_axes = other.m_axes; 
	m_colors = other.m_colors;
	m_centers_calculated = false;
	m_grain_voxel_surfaces_calculated = false;
	m_bad_grain_count=std::string("");
	m_centers.clear();

	calculate_surfaces_loop_spawn(config);
	m_data_voxel.simple_volume_finder(other.m_data_voxel.get_a_voxel_const_ref(),config,other.m_data_voxel.how_many_grains().first);
}

void VoxelEnvelope::grow_them_a_little(VoxelEnvelope& other,VoxelEnvelope& source,Config& config,std::string name)
{
	m_myname = name;
	m_axes = other.m_axes; 
	m_colors = other.m_colors;
	m_centers_calculated = false;
	m_grain_voxel_surfaces_calculated = false;
	m_bad_grain_count=std::string("");
	m_centers.clear();

	calculate_surfaces_loop_spawn(config);
	m_data_voxel.simple_volume_grow(other.m_data_voxel.get_a_voxel_const_ref(),source.m_data_voxel.get_a_voxel_const_ref());
}
		
void VoxelEnvelope::calculate_difference(VoxelEnvelope& a,VoxelEnvelope& b,Config& config,std::string name)
{
	m_myname = name;
	m_axes = a.m_axes; 
	m_colors = a.m_colors;
	m_centers_calculated = false;
	m_grain_voxel_surfaces_calculated = false;
	m_bad_grain_count=std::string("");
	m_centers=a.m_centers;

	m_data_voxel.diff(a.m_data_voxel.get_a_voxel_const_ref(),b.m_data_voxel.get_a_voxel_const_ref());
	calculate_grain_voxel_surfaces(config);
}

void VoxelEnvelope::calculate_final_version(VoxelEnvelope& a,VoxelEnvelope& b,Config& config,std::string name)
{
	m_myname = name;
	m_axes = a.m_axes; 
	m_colors = a.m_colors;
	m_centers_calculated = false;
	m_grain_voxel_surfaces_calculated = false;
	m_bad_grain_count=std::string("");
	m_centers=a.m_centers;

	m_data_voxel.final(a.m_data_voxel.get_a_voxel_const_ref(),b.m_data_voxel.get_a_voxel_const_ref(),config);
	calculate_grain_voxel_surfaces(config);
}

void VoxelEnvelope::cut_off_all_tentacles(VoxelEnvelope& other,std::string name,Config& config)
{
	m_myname = name;
	m_axes = other.m_axes; 
	m_colors = other.m_colors;
	m_centers_calculated = false;
	m_grain_voxel_surfaces_calculated = false;
	m_bad_grain_count=std::string("");
	m_centers=other.m_centers;

//	calculate_grain_voxel_surfaces(config);
//	m_centers=other.m_centers;
	m_data_voxel.cut_off_junk(other.m_data_voxel.get_a_voxel_const_ref(),other.m_data_voxel.how_many_grains().first,m_centers);
	calculate_grain_voxel_surfaces(config);
//	calculate_surfaces_loop_spawn(config);
}

void VoxelEnvelope::draw(Config& c)
{
	m_data_voxel.draw(m_axes,m_colors,c);
	m_data_surface.draw(m_axes,m_colors,c,m_centers);
}

void VoxelEnvelope::calculate_centers(Config& c)
{
	if(centers_calculated())
		return;
	bool ready_=ready();
	t_voxel_data dat(m_data_voxel.get_a_voxel_copy());

	int si=dat.size();
	int sx=dat[0].size();
	int sy=dat[0][0].size();

	std::vector<Vector3r> centers;
	centers.resize(c.grains()+1,Vector3r(0,0,0));
	std::vector<float> counts;
	counts.resize(c.grains()+1,0.0);

	for(int i=0;i<si;i++)
		for(int x=0;x<sx;x++)
			for(int y=0;y<sy;y++)
			{
				int gr=dat[i][x][y];
				if(gr!=0)
				{
					centers[gr]+=Vector3r(x,y,i);
					counts[gr]+=1.0;
				}
			}

	for(int i=0;i<c.grains()+1;++i)
		if(counts[i]>0)
		{
			centers[i]/=counts[i];
			std::cout << "Center of grain " << i << ": " << centers[i] << "\n";
		}

	boost::mutex::scoped_lock scoped_lock(m_voxel_envelope_mutex);
	if(m_centers_calculated)
		return;
	m_centers=centers;
	if(ready_) // so that was the final version
		m_centers_calculated=true;
}

void VoxelEnvelope::calculate_grain_voxel_surfaces(Config& c)
{
	if(grain_voxel_surfaces_calculated())
		return;
	calculate_centers(c);
	bool ready_ = ready();
	bool centers_calculated_ = centers_calculated();

	std::pair<std::set<int>,std::string> tmp = m_data_voxel.how_many_grains();
	int how_many_grains_there_are = tmp.first.size();
	if(how_many_grains_there_are != c.grains())
		m_bad_grain_count = " \nThere are "+boost::lexical_cast<std::string>(how_many_grains_there_are)
		+ " grains in the sample, not "+boost::lexical_cast<std::string>(c.grains())+" as specified!!\n"
		+ "Missing grains: " + tmp.second + "\n";

	DataSurface data_surface;
	data_surface.create_surfaces(m_data_voxel.get_a_voxel_copy() , tmp.first , this );

	boost::mutex::scoped_lock scoped_lock(m_voxel_envelope_mutex);
	if(m_grain_voxel_surfaces_calculated)
		return;
	m_data_surface = data_surface;
	if(ready_ && centers_calculated_)
		m_grain_voxel_surfaces_calculated=true;
}

Vector3r VoxelEnvelope::center(Config& c)
{
	Vector3r res(0,0,0);
	calculate_centers(c);
	float size=c.selection_copy().size();
	if(size>0)
	{
		BOOST_FOREACH(int i,c.selection_copy())
			res+=m_centers[i];
		res/=size;
	}
	return res;
}
		
std::string VoxelEnvelope::message()
{
	std::string res;
	typedef std::string s;
	res+=(ready()?s("data loaded"):s("still loading data"))+".\n";
	res+="grain centers "+(centers_calculated()?s(""):s("not "))+"calculated.\n";
	res+="grain voxel surfaces "+(grain_voxel_surfaces_calculated()?s(""):s("not "))+"calculated.\n";
	res+=m_bad_grain_count;
	return res;
}

bool VoxelEnvelope::centers_calculated()
{
	boost::mutex::scoped_lock scoped_lock(m_voxel_envelope_mutex);
	return m_centers_calculated;
}

bool VoxelEnvelope::grain_voxel_surfaces_calculated()
{
	boost::mutex::scoped_lock scoped_lock(m_voxel_envelope_mutex);
	return m_grain_voxel_surfaces_calculated;
}

bool VoxelEnvelope::ready()
{
	return m_data_voxel.ready();
}


YADE_REQUIRE_FEATURE(PHYSPAR);

