YADE_REQUIRE_FEATURE(COMPILED_BY_INCLUDE_FROM_SnowVoxelsLoader.cpp)
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "DataVoxel.hpp"
#include "Config.hpp"
#include <list>
#include <algorithm>
#include <vector>
#include <stdlib.h>

DataVoxel::DataVoxel() : m_ready(false)
{
	m_voxel_data.clear();
}

void DataVoxel::draw(SafeVectors3& axes,SafeVectors3& colors,Config& c)
{
	if(c.draw_using_surfaces())
		return;
	boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	int DD=c.display_steps();
	bool selected=c.selected();
	bool names=c.draw_with_names();
	boost::mutex::scoped_lock lock(m_voxel_mutex);
	int si=m_voxel_data.size();
	if(si==0)
		return;
	int sx=m_voxel_data[0].size();
	int sy=m_voxel_data[0][0].size();
	if(!names)
		c.apply_transparency();
	for(int i=0;i<si;i+=DD)
	{
		for(int x=0;x<sx;x+=DD)
			for(int y=0;y<sy;y+=DD)
			{
				int gr=m_voxel_data[i][x][y];
				if( (selected && c.display_selection()) ? c.selected(gr) : gr!=0 )
				{
					if(!names)
					{
						Vector3r co=colors[gr];//-Vector3r(1,1,1)*0.07;
						if(c.transparency())
							glColor4f(co[0],co[1],co[2],0.02);
						else
							glColor3f(co[0],co[1],co[2]);
						if(!c.display_selection())
							c.light_selection(gr);
					}
					if(names) glPushName(gr);
					glTranslatef(x,y,i);
					glutSolidCube(DD);
					glTranslatef(-x,-y,-i);
					if(names) glPopName();
				}
			}
		if(c.time_draw_limit() && 
		   boost::posix_time::time_duration(boost::posix_time::microsec_clock::local_time() - start) > boost::posix_time::milliseconds(400))
			return;
	}
	if(!names)
		c.remove_transparency();
}

void DataVoxel::save_txt_result(Config& c)
{
	for(int i=c.start() ; i<c.end()+1 ; ++i)
	{
		std::string number(boost::lexical_cast<std::string>(i));
		while(number.length()<c.digits()) number="0"+number;
		std::string fname=c.file_basename()+number+".txt";
	
		std::ofstream out(fname.c_str());
	
		std::cerr << "saving " << fname << "... ";
		
		int sx = m_voxel_data[0].size();
		if(sx == 0)
			std::cerr << "sx==0, wtf?\n";
		int sy = m_voxel_data[0][0].size();
		if(sy == 0)
			std::cerr << "sy==0, wtf?\n";
		out << (int)sx << " " << (int)sy << " \n";
		for(int y=0;y<sy;++y)
		{
			for(int x=0;x<sx;++x)
			{
				out << (int)m_voxel_data[c.end()-i][y][x] << " ";
			}
			out << "\n";
		}
		out.close();
		std::cerr << "done.\n";
	}
}

void DataVoxel::read_txt(std::string fname,t_voxel_slice& dat)
{
	std::cerr << "loading " << fname << "... ";
	std::ifstream file(fname.c_str());
	if(!file)
		std::cerr << "\nERROR: file " << fname <<" could not be opened for reading\n", exit(0);

	int sx,sy,n;
	file >> sx >> sy;
	dat.resize(sy);
	for(int j=0;j<sy;++j)
	{
		dat[j].resize(sx,0);
		for(int i=0;i<sx;++i)
		{
			file >> n;

			if(n > 254)
				std::cerr << "WARNING: To save memory only 254 grains is allowed! Modify typedef t_voxel in DataVoxel.hpp:8.\n";

			dat[j][i]=n;
		}
	}
	file.close();
	std::cerr << "loaded.\n";
}

void DataVoxel::load(Config& c)
{
	if(ready())
		return;
	//m_voxel_data.clear();
	for(int i=c.start() ; i<c.end()+1 ; ++i)
	{
		t_voxel_slice slice;
		std::string number(boost::lexical_cast<std::string>(i));
		while(number.length()<c.digits()) number="0"+number;
		read_txt(std::string(c.file_basename()+number+".txt"),slice);

		{
			boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
			m_voxel_data.insert(m_voxel_data.begin(),slice);
		}
	}

	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	m_ready=true;
}

int  DataVoxel::volume_of_grain_id(const t_voxel_data& data, int GRAIN_ID)
{
	int result=0;
	for(size_t i = 0 ; i < data.size() ; ++i)
		for( size_t x = 0 ; x < data[i].size() ; ++x)
			for( size_t y = 0 ; y < data[i][x].size() ; ++y)
				if(data[i][x][y] == GRAIN_ID)
					++result;
	return result;
}

void DataVoxel::simple_volume_grow(const t_voxel_data& other,const t_voxel_data& source)
{
	if(ready())
		return;

	{
		boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
		m_voxel_data = source;
	}
	int d_height = other.size();
	int d_width = other[0].size();
	int d_depth = other[0][0].size();
	// now add volumes separated from other grains
	bool anything=true;
	while(anything)
	{
		anything=false;
		for(int j = 0; j < d_height; j++)
		for(int i = 0; i < d_width; i++)
		for(int h = 0; h < d_depth; h++)
			if((m_voxel_data[j][i][h] == 0) && (other[j][i][h] != 0))
			{
if((j-1>0       ) && (m_voxel_data[j-1][i  ][h  ] != 0)) {boost::mutex::scoped_lock scoped_lock(m_voxel_mutex); m_voxel_data[j][i][h] = m_voxel_data[j-1][i  ][h  ]; anything=true;} else
if((i-1>0       ) && (m_voxel_data[j  ][i-1][h  ] != 0)) {boost::mutex::scoped_lock scoped_lock(m_voxel_mutex); m_voxel_data[j][i][h] = m_voxel_data[j  ][i-1][h  ]; anything=true;} else
if((i+1<d_width ) && (m_voxel_data[j  ][i+1][h  ] != 0)) {boost::mutex::scoped_lock scoped_lock(m_voxel_mutex); m_voxel_data[j][i][h] = m_voxel_data[j  ][i+1][h  ]; anything=true;} else
if((j+1<d_height) && (m_voxel_data[j+1][i  ][h  ] != 0)) {boost::mutex::scoped_lock scoped_lock(m_voxel_mutex); m_voxel_data[j][i][h] = m_voxel_data[j+1][i  ][h  ]; anything=true;} else
if((h+1<d_depth ) && (m_voxel_data[j  ][i  ][h+1] != 0)) {boost::mutex::scoped_lock scoped_lock(m_voxel_mutex); m_voxel_data[j][i][h] = m_voxel_data[j  ][i  ][h+1]; anything=true;} else
if((h-1>0       ) && (m_voxel_data[j  ][i  ][h-1] != 0)) {boost::mutex::scoped_lock scoped_lock(m_voxel_mutex); m_voxel_data[j][i][h] = m_voxel_data[j  ][i  ][h-1]; anything=true;};
			}
	}
	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	m_ready=true;
};

void DataVoxel::diff(const t_voxel_data& a,const t_voxel_data& b)
{
	if(ready())
		return;

	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	m_voxel_data = a;

	int d_height = a.size();
	int d_width = a[0].size();
	int d_depth = a[0][0].size();
	
	for(int j = 0; j < d_height; j++)
	for(int i = 0; i < d_width; i++)
	for(int h = 0; h < d_depth; h++)
	{
		m_voxel_data[j][i][h] = 0;

//		if(((a[j][i][h] == 0) && (b[j][i][h] != 0))) 
//			m_voxel_data[j][i][h] = b[j][i][h];
//		if(((a[j][i][h] != 0) && (b[j][i][h] == 0)))
//			m_voxel_data[j][i][h] = a[j][i][h];

		if( a[j][i][h] != b[j][i][h] )
			m_voxel_data[j][i][h] = (a[j][i][h] == 0 ? 40 : a[j][i][h] );
	}
	m_ready=true;
};

		
void DataVoxel::final(const t_voxel_data& experimental,const t_voxel_data& clean,Config& config)
{
	std::vector<Vector3<int> > needs_to_be_checked;

	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	m_voxel_data = clean;

	for(int id = 1 ; id<config.grains()+10 ; ++id)
	{	
//		std::cerr << "checking grain " << id << " for clarity\n";
		if( volume_of_grain_id(clean,id) < volume_of_grain_id(experimental,id)/5 )
		{
			// clean this one off
			std::cerr << "grain " << id << " too small - removing\n";
			for(size_t i = 0 ; i < m_voxel_data.size() ; ++i)
				for( size_t x = 0 ; x < m_voxel_data[i].size() ; ++x)
					for( size_t y = 0 ; y < m_voxel_data[i][x].size() ; ++y)
						if(m_voxel_data[i][x][y] == id)
							m_voxel_data[i][x][y] = 0;
		}
	}

	needs_to_be_checked.clear();
	for(size_t i = 0 ; i < m_voxel_data.size() ; ++i)
		for( size_t x = 0 ; x < m_voxel_data[i].size() ; ++x)
			for( size_t y = 0 ; y < m_voxel_data[i][x].size() ; ++y)
				if(m_voxel_data[i][x][y] != experimental[i][x][y] && m_voxel_data[i][x][y] == 0)
					needs_to_be_checked.push_back(Vector3<int>(i,x,y));

	std::random_shuffle(needs_to_be_checked.begin(),needs_to_be_checked.end());
	std::vector<int> checked;
	checked.resize(needs_to_be_checked.size(),0);
	int checked_end = 0;
	int last_checked = 0;
	int dead = 0;
	bool anything=true;
	std::vector<int>  s;
	s.resize(config.grains()+10,0);
	while(anything && dead < 10)
	{
		anything=false;
		BOOST_FOREACH(int c,checked)
		{
			if(c == 0)
			{
				anything=true;
				std::cerr << "finalizing (dead: "<<dead<<"): " << last_checked << "," << checked_end << "/" << checked.size() << "\n";
				break;
			}
		}
		last_checked = checked_end;
		for(size_t D = 0 ; D < needs_to_be_checked.size() ; ++D)
		{
			if(checked[D] != 0)
				continue;

			int i = needs_to_be_checked[D][0];
			int x = needs_to_be_checked[D][1];
			int y = needs_to_be_checked[D][2];
			if( 
				m_voxel_data[i+1][x  ][y  ] != 0 ||
				m_voxel_data[i-1][x  ][y  ] != 0 ||
				m_voxel_data[i  ][x+1][y  ] != 0 ||
				m_voxel_data[i  ][x-1][y  ] != 0 ||
				m_voxel_data[i  ][x  ][y+1] != 0 ||
				m_voxel_data[i  ][x  ][y-1] != 0
			)
			{
				for(size_t ZZ = 0 ; ZZ < s.size() ; ++ZZ) 
					s[ZZ] = 0;

				if(m_voxel_data[i+1][x  ][y  ] != 0){ s[m_voxel_data[i+1][x  ][y  ]] += 1; };
				if(m_voxel_data[i-1][x  ][y  ] != 0){ s[m_voxel_data[i-1][x  ][y  ]] += 1; };
				if(m_voxel_data[i  ][x+1][y  ] != 0){ s[m_voxel_data[i  ][x+1][y  ]] += 1; };
				if(m_voxel_data[i  ][x-1][y  ] != 0){ s[m_voxel_data[i  ][x-1][y  ]] += 1; };
				if(m_voxel_data[i  ][x  ][y+1] != 0){ s[m_voxel_data[i  ][x  ][y+1]] += 1; };
				if(m_voxel_data[i  ][x  ][y-1] != 0){ s[m_voxel_data[i  ][x  ][y-1]] += 1; };

				int max_id=0;
				int max_contact = 0;
				for(size_t ZZ = 0 ; ZZ < s.size() ; ++ZZ)
				{
					if(s[ZZ] > max_contact)
					{
						max_id = ZZ;
						max_contact = s[max_id];
					}
					if(max_id != 0 && s[ZZ] == max_contact && rand()%2 == 0)
					{
						max_id = ZZ;
						max_contact = s[max_id];
					}
				}

				if(max_id == 0)
					std::cerr << "WTF - max_id==0 ??\n";
				if(m_voxel_data[i][x][y] != 0)
					std::cerr << "WTF nonzero ??\n";
				m_voxel_data[i][x][y] = max_id;
				checked[D] = 1;
				++checked_end;
			}

/*			{
				m_voxel_data[i][x][y] = 50;
				checked[D] = 1;
				++checked_end;
			}
*/		}
		if(last_checked == checked_end)
			++dead;
	}

	std::set<int> removed;removed.clear();
	for(int id = 1 ; id<config.grains()+10 ; ++id)
	{// change owners due to overhelming contact surface
		std::vector<int>  surfaces;
		surfaces.resize(config.grains()+10,0);
		for(size_t ZZ = 0 ; ZZ < surfaces.size() ; ++ZZ) 
			surfaces[ZZ] = 0;
		int total_surface = 0;

		for(size_t i = 0 ; i < m_voxel_data.size() ; ++i)
			for( size_t x = 0 ; x < m_voxel_data[i].size() ; ++x)
				for( size_t y = 0 ; y < m_voxel_data[i][x].size() ; ++y)
					if(m_voxel_data[i][x][y] == id)
					{
						// calculate surface of this grain, and surface of contact with other grains
						for(int I = -1 ; I < 2 ; ++I)
							for(int X = -1 ; X < 2 ; ++X)
								for(int Y = -1 ; Y < 2 ; ++Y)
									if(   i+I>=0 && i+I<m_voxel_data.size()
									   && x+X>=0 && x+X<m_voxel_data[i].size()
									   && y+Y>=0 && y+Y<m_voxel_data[i][x].size())
									{
										surfaces[m_voxel_data[i+I][x+X][y+Y]] += 1;

										if( m_voxel_data[i+I][x+X][y+Y] != id )
											total_surface += 1;
									}
					}
		int max_id=0;
		int max_contact = 0;
		for(size_t ZZ = 1 ; ZZ < surfaces.size() ; ++ZZ)
		{
			if(surfaces[ZZ] > max_contact && ZZ != (unsigned int)(id))
			{
				max_id = ZZ;
				max_contact = surfaces[max_id];
			}
		}
		if( max_contact > total_surface*0.3 && removed.find(id) == removed.end())
		{
			std::cerr << "grain " << id << " has contact with grain "<< max_id <<" bigger than 0.3 of its total surface ( " << max_contact << "/" << total_surface << " = "<<100*max_contact/total_surface<<" \% ) - removing\n";
			// clean this one off
			removed.insert(max_id);
			for(size_t i = 0 ; i < m_voxel_data.size() ; ++i)
				for( size_t x = 0 ; x < m_voxel_data[i].size() ; ++x)
					for( size_t y = 0 ; y < m_voxel_data[i][x].size() ; ++y)
						if(m_voxel_data[i][x][y] == id)
							m_voxel_data[i][x][y] = max_id;
		}
		else if(total_surface > 0)
			std::cerr << "grain " << id << " has contact with grain "<< max_id <<" percentage " << 100*max_contact/total_surface << " \% of its total surface (" << max_contact << "/" << total_surface << ")\n";
		else	std::cerr << "grain " << id << " has zero total_surface!\n";
	}
	m_ready=true;
}

void DataVoxel::cut_off_junk(const t_voxel_data& other,std::set<int> grains,const std::vector<Vector3r>& centers)
{
	if(ready())
		return;

	{
		boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
		// first resize to other's size
		m_voxel_data.resize(other.size());
		for(size_t i = 0 ; i < m_voxel_data.size() ; ++i)
		{
			m_voxel_data[i].resize(other[0].size());
			for( size_t x = 0 ; x < m_voxel_data[i].size() ; ++x)
			{
				m_voxel_data[i][x].resize(other[0][0].size() , 0);
				for( size_t y = 0 ; y < m_voxel_data[i][x].size() ; ++y)
					m_voxel_data[i][x][y] = 0;
			}
		}
	}

	std::vector<std::vector<std::vector<unsigned char> > > box;
	int d_height = other.size();
	int d_width  = other[0].size();
	int d_depth  = other[0][0].size();
	box.resize(d_height);
	for(int j = 0; j < d_height; j++)
	{
		box[j].resize(d_width);
		for(int i = 0 ; i <d_width ; ++i)
		{
			box[j][i].resize(d_depth,0);
			for(int h = 0; h < d_depth; h++)
			{
				box[j][i][h] = 0;
			}
		}
	}
	BOOST_FOREACH(int grain_id, grains)
	{
		if(grain_id != 0)
		{
			std::cerr << "cut off tentacles, grain: " << grain_id << "\n";

			//std::cerr << centers[grain_id] << "\n";
			if(other[centers[grain_id][2]][centers[grain_id][0]][centers[grain_id][1]] == grain_id)
			{
				grow_balloon(centers[grain_id],grain_id,other,box);
			}
			else
			{
				std::cerr << "grain needs a better center: " << grain_id << "\n";
				int J = centers[grain_id][2];
				int I = centers[grain_id][0];
				int K = centers[grain_id][1];
				Vector3r cen(0,0,0);

				int i_up=0,i_dn=0,j_up=0,j_dn=0,k_up=0,k_dn=0;
				while((J+j_up <  d_height-1) && (other[J+j_up][I     ][K     ] != grain_id)) ++j_up;
				while((J-j_dn >  0         ) && (other[J-j_dn][I     ][K     ] != grain_id)) ++j_dn;
				while((I+i_up <  d_width-1 ) && (other[J     ][I+i_up][K     ] != grain_id)) ++i_up;
				while((I-i_dn >  0         ) && (other[J     ][I-i_dn][K     ] != grain_id)) ++i_dn;
				while((K+k_up <  d_depth -1) && (other[J     ][I     ][K+k_up] != grain_id)) ++k_up;
				while((K-k_dn >  0         ) && (other[J     ][I     ][K-k_dn] != grain_id)) ++k_dn;

				int shortest = std::min( std::min(i_up,i_dn) , std::min( std::min(j_up,j_dn) , std::min(k_up,k_dn) ) );

				if(shortest == j_up && J+j_up <  d_height)
					cen =			Vector3r( I      , K      , J+j_up );
				else
				if(shortest == j_dn && J-j_dn >= 0       )
					cen =			Vector3r( I      , K      , J+j_dn );
				else
				if(shortest == i_up && I+i_up <  d_width )
					cen =			Vector3r( I+i_up , K      , J      );
				else
				if(shortest == i_dn && I-i_dn >= 0       )
					cen =			Vector3r( I-i_dn , K      , J      );
				else
				if(shortest == k_up && K+k_up <  d_depth )
					cen =			Vector3r( I      , K+k_up , J      );
				else
				if(shortest == k_dn && K-k_dn >= 0       )
					cen =			Vector3r( I      , K-k_dn , J      );
				else
				{
					std::cerr << "Where is the center?????\n";
				}
				if(other[cen[2]][cen[0]][cen[1]] == grain_id)
					grow_balloon(cen,grain_id,other,box);
				else
					std::cerr << "Totally wtf.\n";
			}
		}
	}
	std::cerr << "cutting tentacles finished.\n";
	
	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	m_ready=true;
};

void DataVoxel::grow_balloon(Vector3r center, int GRAIN_ID,const t_voxel_data& other, std::vector<std::vector<std::vector<unsigned char> > >& box)
{
//	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);

	int d_height = other.size();
	int d_width  = other[0].size();
	int d_depth  = other[0][0].size();

	int min_j = 10000, max_j = -1;
	int min_i = 10000, max_i = -1;
	int min_h = 10000, max_h = -1;

	for(int j = 0; j < d_height; j++)
		for(int i = 0 ; i <d_width ; ++i)
			for(int h = 0; h < d_depth; h++)
			{
				if(other[j][i][h] == GRAIN_ID)
				{
					min_j = std::min(min_j,j); max_j = std::max(max_j,j);
					min_i = std::min(min_i,i); max_i = std::max(max_i,i);
					min_h = std::min(min_h,h); max_h = std::max(max_h,h);
				}
			}

	m_voxel_data[center[2]][center[0]][center[1]] = GRAIN_ID;
	std::list<Vector3<int> > cur,nowe;
	nowe.clear();
	nowe.push_back(Vector3<int>(center[2],center[0],center[1]));

	bool anything=true;
	while(anything)
	{
		cur = nowe;
		nowe.clear();

		anything=false;
		BOOST_FOREACH(Vector3<int> vvv, cur)
		{
			int cj = vvv[0];
			int ci = vvv[1];
			int ch = vvv[2];
			for(int zj = -1 ; zj <= 1 ; ++zj)
				for(int zi = -1 ; zi <= 1 ; ++zi)
					for(int zh = -1 ; zh <= 1 ; ++zh)
					{
						int j = cj+zj;
						int i = ci+zi;
						int h = ch+zh;
						if((!(zj==0 && zi==0 && zh == 0)) &&
						   (j-1>0       ) && 
						   (i-1>0       ) && 
						   (i+1<d_width ) && 
						   (j+1<d_height) && 
						   (h+1<d_depth ) && 
						   (h-1>0       )
						  )
						{
							if((m_voxel_data[j][i][h] == 0) && (other[j][i][h] == GRAIN_ID) && (box[j][i][h] == 0) &&
							   (
							    ((m_voxel_data[j-1][i  ][h  ] == GRAIN_ID)) ||
							    ((m_voxel_data[j  ][i-1][h  ] == GRAIN_ID)) ||
							    ((m_voxel_data[j  ][i+1][h  ] == GRAIN_ID)) ||
							    ((m_voxel_data[j+1][i  ][h  ] == GRAIN_ID)) ||
							    ((m_voxel_data[j  ][i  ][h+1] == GRAIN_ID)) ||
							    ((m_voxel_data[j  ][i  ][h-1] == GRAIN_ID))
							   )
							  )
							{
								int NN =0;
								for(int ZJ = -1 ; ZJ <= 1 ; ++ZJ)
									for(int ZI = -1 ; ZI <= 1 ; ++ZI)
										for(int ZH = -1 ; ZH <= 1 ; ++ZH)
											if((j+ZJ > 0) && (j+ZJ < d_height) &&
											   (i+ZI > 0) && (i+ZI < d_width) &&
											   (h+ZH > 0) && (h+ZH < d_depth) &&
											   (other[j+ZJ][i+ZI][h+ZH] == GRAIN_ID))
											++NN;
								if(NN == 9)
								{
									boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
									m_voxel_data[j][i][h] = GRAIN_ID;
									anything = true;
									box[j][i][h] = 1;
									nowe.push_back(Vector3<int>(j,i,h));
									//std::cerr << "fast 9\n";
								}
								else
								{
									int SIZE = 23;//35;//15;//31;
									int RANGE = (SIZE-1)/2;
									int VOL = SIZE*SIZE*SIZE;
									int neighbours = 0;
									Vector3r vec(0,0,0);
									float count = 0;
									for(int J = -1*RANGE ; J <= RANGE ; ++J)
										for(int I = -1*RANGE ; I <= RANGE ; ++I)
											for(int H = -1*RANGE ; H <= RANGE ; ++H)
												if((j+J > 0) && (j+J < d_height) &&
												   (i+I > 0) && (i+I < d_width ) &&
												   (h+H > 0) && (h+H < d_depth ) &&
												   (other[j+J][i+I][h+H] == GRAIN_ID))
												{
													++neighbours;
													vec += Vector3r(J,I,H);
													count +=1.0;
												}
									if(count>=1.0)
									{
										float len = vec.Length()/count;
										if(neighbours >= VOL*4/10 && len > RANGE*0.40)
										{
											boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
											m_voxel_data[j][i][h] = GRAIN_ID;
											anything = true;
											box[j][i][h] = 1;
											nowe.push_back(Vector3<int>(j,i,h));
										} else
										if(neighbours >= VOL*5/10 && len > RANGE*0.25)
										{
											boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
											m_voxel_data[j][i][h] = GRAIN_ID;
											anything = true;
											box[j][i][h] = 1;
											nowe.push_back(Vector3<int>(j,i,h));
										} else
										if(neighbours >= VOL*7/10 && len > RANGE*0.05)
										{
											boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
											m_voxel_data[j][i][h] = GRAIN_ID;
											anything = true;
											box[j][i][h] = 1;
											nowe.push_back(Vector3<int>(j,i,h));
										} else
										if(neighbours >= VOL*9/10)
										{
											boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
											m_voxel_data[j][i][h] = GRAIN_ID;
											anything = true;
											box[j][i][h] = 1;
											nowe.push_back(Vector3<int>(j,i,h));
										} else
											box[j][i][h] = 2;
									} else
										box[j][i][h] = 4;
								}
							}
							else
							{
								if(
								   ((j-1>0       ) && (m_voxel_data[j-1][i  ][h  ] == GRAIN_ID)) ||
								   ((i-1>0       ) && (m_voxel_data[j  ][i-1][h  ] == GRAIN_ID)) ||
								   ((i+1<d_width ) && (m_voxel_data[j  ][i+1][h  ] == GRAIN_ID)) ||
								   ((j+1<d_height) && (m_voxel_data[j+1][i  ][h  ] == GRAIN_ID)) ||
								   ((h+1<d_depth ) && (m_voxel_data[j  ][i  ][h+1] == GRAIN_ID)) ||
								   ((h-1>0       ) && (m_voxel_data[j  ][i  ][h-1] == GRAIN_ID))
								  )
								{
									box[j][i][h] = 3;
								}
							}
						}
					}
		}
	}
}

void DataVoxel::simple_volume_finder(const t_voxel_data& other,Config& config,std::set<int> grains)
{
	if(ready())
		return;

	int d_height = other.size();
	int d_width = other[0].size();
	int d_depth = other[0][0].size();
	{
		boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
		// first resize to other's size
		m_voxel_data.resize(other.size());
		for(size_t i = 0 ; i < m_voxel_data.size() ; ++i)
		{
			m_voxel_data[i].resize(other[0].size());
			for( size_t x = 0 ; x < m_voxel_data[i].size() ; ++x)
			{
				m_voxel_data[i][x].resize(other[0][0].size() , 0);
				for( size_t y = 0 ; y < m_voxel_data[i][x].size() ; ++y)
					m_voxel_data[i][x][y] = 0;
			}
		}
	}

	std::vector<std::vector<std::vector<unsigned char> > > box;
	BOOST_FOREACH(int grain_id, grains)
	{
		int how_many=0;
		int biggest=-1;
		if(grain_id != 0)
		{ // we don't want grain 0 - it's empty space around grains.
			std::cerr << "simple volume finder, grain: " << grain_id << "\n";

			find_volumes(other, box, how_many, grain_id);
			// now box has all separate volumes of selected grain_id
			if(how_many != 1)
			{ // we need to pick the biggest one
				std::vector<int> counts;
				counts.resize(how_many+2,0);
				for(int j = 0; j < d_height; j++)
					for(int i = 0; i < d_width; i++)
						for(int h = 0; h < d_depth; h++)
							if(box[j][i][h] != 0)
								counts[box[j][i][h]] += 1;

				// determine the biggest volume
				int max_count = 0;
				for(int h = 0; h < how_many+2; ++h)
				{
					if(counts[h] > max_count)
					{
						max_count = counts[h];
						biggest = h;
					}
				}
			}
			else if(how_many == 1)
			{
				biggest = 1;
			}

			if(biggest != -1)
			{ // ok, we got the biggest volume, we need to "copy" it into our local m_voxel_data
				boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
				for(int j = 0; j < d_height; j++)
					for(int i = 0; i < d_width; i++)
						for(int h = 0; h < d_depth; h++)
							if(box[j][i][h] == biggest)
								m_voxel_data[j][i][h] = grain_id;
			}
			else
			{
				std::cerr << "\nCannot find biggest volume, WTF?? "<< grain_id <<".\n";
			}
		}
	}
	std::cerr << "finding grain volumes finished.\n";

	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	m_ready=true;
};


void DataVoxel::find_volumes(const t_voxel_data& data, std::vector<std::vector<std::vector<unsigned char> > >& box, int& how_many, int GRAIN_ID)
{
	//void MapGenerator::findAreasOf(Tile::Type THIS_TILE,std::vector<std::vector<int> >& box,int& how_many)

	int d_height = data.size();
	int d_width = data[0].size();
	int d_depth = data[0][0].size();

	how_many = 0;

	int min_j = 10000, max_j = -1;
	int min_i = 10000, max_i = -1;
	int min_h = 10000, max_h = -1;

	box.resize(d_height);
	for(int j = 0; j < d_height; j++)
	{
		box[j].resize(d_width);
		for(int i = 0 ; i <d_width ; ++i)
		{
			box[j][i].resize(d_depth,0);
			for(int h = 0; h < d_depth; h++)
			{
				box[j][i][h] = 0;
				if(data[j][i][h] == GRAIN_ID)
				{
					min_j = std::min(min_j,j); max_j = std::max(max_j,j);
					min_i = std::min(min_i,i); max_i = std::max(max_i,i);
					min_h = std::min(min_h,h); max_h = std::max(max_h,h);
				}
			}
		}
	}
	min_j = std::max(min_j , 1); max_j = std::min(max_j , d_height-1);
	min_i = std::max(min_i , 1); max_i = std::min(max_i , d_width-1 );
	min_h = std::max(min_h , 1); max_h = std::min(max_h , d_depth-1 );

	std::cerr 
	<< min_j << " " << max_j << " " << d_height<< "\n"
	<< min_i << " " << max_i << " " << d_width << "\n" 
	<< min_h << " " << max_h << " " << d_depth << "\n"; 


	// find all enclosed areas by scanning the map
	// distinct areas have different numbers in box
	for(int j = min_j; j < max_j; j++)
		for(int i = min_i; i < max_i; i++)
			for(int h = min_h; h < max_h; h++)
				if (box[j][i][h]==0 &&
						data[j][i][h] == GRAIN_ID &&
						(
						 (
						  data[(j-1)][i-1][h  ] == GRAIN_ID &&
						  data[(j  )][i-1][h  ] == GRAIN_ID &&
						  data[(j-1)][i  ][h  ] == GRAIN_ID &&

						  data[(j-1)][i-1][h+1] == GRAIN_ID &&
						  data[(j  )][i-1][h+1] == GRAIN_ID &&
						  data[(j-1)][i  ][h+1] == GRAIN_ID &&
						  data[(j  )][i  ][h+1] == GRAIN_ID
						 )  ||

						 (
						  data[(j-1)][i-1][h  ] == GRAIN_ID &&
						  data[(j  )][i-1][h  ] == GRAIN_ID &&
						  data[(j-1)][i  ][h  ] == GRAIN_ID &&

						  data[(j-1)][i-1][h-1] == GRAIN_ID &&
						  data[(j  )][i-1][h-1] == GRAIN_ID &&
						  data[(j-1)][i  ][h-1] == GRAIN_ID &&
						  data[(j  )][i  ][h-1] == GRAIN_ID
						 )  ||

						 (
						  data[(j-1)][i  ][h  ] == GRAIN_ID &&
						  data[(j-1)][i+1][h  ] == GRAIN_ID &&
						  data[(j  )][i+1][h  ] == GRAIN_ID &&

						  data[(j-1)][i  ][h+1] == GRAIN_ID &&
						  data[(j-1)][i+1][h+1] == GRAIN_ID &&
						  data[(j  )][i+1][h+1] == GRAIN_ID &&
						  data[(j  )][i  ][h+1] == GRAIN_ID
						 ) ||

						 (
						  data[(j-1)][i  ][h  ] == GRAIN_ID &&
						  data[(j-1)][i+1][h  ] == GRAIN_ID &&
						  data[(j  )][i+1][h  ] == GRAIN_ID &&

						  data[(j-1)][i  ][h-1] == GRAIN_ID &&
						  data[(j-1)][i+1][h-1] == GRAIN_ID &&
						  data[(j  )][i+1][h-1] == GRAIN_ID &&
						  data[(j  )][i  ][h-1] == GRAIN_ID
						 ) ||

						 (
						  data[(j  )][i+1][h  ] == GRAIN_ID &&
						  data[(j+1)][i+1][h  ] == GRAIN_ID &&
						  data[(j+1)][i  ][h  ] == GRAIN_ID &&

						  data[(j  )][i+1][h+1] == GRAIN_ID &&
						  data[(j+1)][i+1][h+1] == GRAIN_ID &&
						  data[(j+1)][i  ][h+1] == GRAIN_ID &&
						  data[(j  )][i  ][h+1] == GRAIN_ID
						 ) ||

						 (
						  data[(j  )][i+1][h  ] == GRAIN_ID &&
						  data[(j+1)][i+1][h  ] == GRAIN_ID &&
						  data[(j+1)][i  ][h  ] == GRAIN_ID &&

						  data[(j  )][i+1][h-1] == GRAIN_ID &&
						  data[(j+1)][i+1][h-1] == GRAIN_ID &&
						  data[(j+1)][i  ][h-1] == GRAIN_ID &&
						  data[(j  )][i  ][h-1] == GRAIN_ID
						 ) ||

						 (
						  data[(j+1)][i  ][h  ] == GRAIN_ID &&
						  data[(j+1)][i-1][h  ] == GRAIN_ID &&
						  data[(j  )][i-1][h  ] == GRAIN_ID &&

						  data[(j+1)][i  ][h+1] == GRAIN_ID &&
						  data[(j+1)][i-1][h+1] == GRAIN_ID &&
						  data[(j  )][i-1][h+1] == GRAIN_ID &&
						  data[(j  )][i  ][h+1] == GRAIN_ID
						 ) ||

						 (
						  data[(j+1)][i  ][h  ] == GRAIN_ID &&
						  data[(j+1)][i-1][h  ] == GRAIN_ID &&
						  data[(j  )][i-1][h  ] == GRAIN_ID &&

						  data[(j+1)][i  ][h-1] == GRAIN_ID &&
						  data[(j+1)][i-1][h-1] == GRAIN_ID &&
						  data[(j  )][i-1][h-1] == GRAIN_ID &&
						  data[(j  )][i  ][h-1] == GRAIN_ID
						 )
						 )
						 )
						 {
							 box[j][i][h]=++how_many+130;
							 int counter=1;
							 while(counter != 0)
							 {
								 counter=0;
								 for(int J = min_j-1; J < max_j+1; J++)
									 for(int I = min_i-1; I < max_i+1; I++)
										 for(int H = min_h-1; H < max_h+1; H++)
										 {
											 if(
												 data[J][I][H] == GRAIN_ID &&
												 box[J][I][H]  == 0 &&
												 (
												  ((J-1>0       ) && (box[J-1][I  ][H  ]==how_many+130)) ||
												  ((I-1>0       ) && (box[J  ][I-1][H  ]==how_many+130)) ||
												  ((I+1<d_width ) && (box[J  ][I+1][H  ]==how_many+130)) ||
												  ((J+1<d_height) && (box[J+1][I  ][H  ]==how_many+130)) ||
												  ((H+1<d_depth ) && (box[J  ][I  ][H+1]==how_many+130)) ||
												  ((H-1>0       ) && (box[J  ][I  ][H-1]==how_many+130))
												 )
											   )
											 {
												 ++counter;
												 box[J][I][H]=how_many+120;
											 }
										 }
								 for(int J = min_j-1; J < max_j+1; J++)
									 for(int I = min_i-1; I < max_i+1; I++)
										 for(int H = min_h-1; H < max_h+1; H++)
										 {
											 if (box[J][I][H]==how_many+130)
												 box[J][I][H]=how_many;
											 if (box[J][I][H]==how_many+120)
												 box[J][I][H]=how_many+130;
										 }
							 }
						 }
	for(int J = 0; J < d_height; J++)
		for(int I = 0; I < d_width; I++)
			for(int H = 0; H < d_depth; H++)
			{
				if (box[J][I][H]>how_many+1)
					std::cerr << "\n DataVoxel::find_volumes - something suspicious here.\n";
			}
};

bool DataVoxel::ready()
{
	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	return m_ready;
};

//void DataVoxel::ready(bool t)
//{
//	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
//	m_ready = t;
//};

t_voxel_data DataVoxel::get_a_voxel_copy()
{
	while(true)
	{
		boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
		if(m_voxel_data.size() > 2)
			return m_voxel_data;
	}
};

std::pair<std::set<int> ,std::string> DataVoxel::how_many_grains()
{
	boost::mutex::scoped_lock scoped_lock(m_voxel_mutex);
	if(m_voxel_data.size() < 2)
	{
		std::cerr << "cannot count grains, must load first\n";
		return std::make_pair(std::set<int>(),"");
	}

	int max=0;
	int si=m_voxel_data.size();
	int sx=m_voxel_data[0].size();
	int sy=m_voxel_data[0][0].size();
	std::set<int> grains;grains.clear();
	for(int i=0;i<si;i++)
	{
		for(int x=0;x<sx;x++)
			for(int y=0;y<sy;y++)
			{
				int gr=m_voxel_data[i][x][y];
				//if( gr!=0 )
				//
				//WITH zero GRAIN too!
					grains.insert(gr);
					max=std::max(max,gr);
			}
	}

	std::string missing("");
	for(int i = 0 ; i<=max ; ++i)
		if(grains.find(i) == grains.end())
			missing += boost::lexical_cast<std::string>(i) + " ";

	return std::make_pair(grains,missing);
};



YADE_REQUIRE_FEATURE(PHYSPAR);

