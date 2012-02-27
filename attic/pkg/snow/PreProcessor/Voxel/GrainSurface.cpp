YADE_REQUIRE_FEATURE(COMPILED_BY_INCLUDE_FROM_SnowVoxelsLoader.cpp)
#include "GrainSurface.hpp"
#include "Config.hpp"

GrainSurface::GrainSurface(t_voxel_data& dat, t_voxel grain_id)
{
	m_min_i = 10000;
	m_min_x = 10000;
	m_min_y = 10000;
	m_max_i = 0;
	m_max_x = 0;
	m_max_y = 0;

	m_own_id = grain_id;
	if(grain_id == 0)
		return;
	
	int si=dat.size();
	int sx=dat[0].size();
	int sy=dat[0][0].size();
	for(t_coord i=0;i<si;++i)
		for(t_coord x=0;x<sx;++x)
			for(t_coord y=0;y<sy;++y)
				if(dat[i][x][y] == grain_id)
				{ // i/x/y
					m_on_i.insert(std::make_pair(x,y)); // x,y  -->  i
					m_on_x.insert(std::make_pair(y,i)); // y,i  -->  x
					m_on_y.insert(std::make_pair(i,x)); // i,x  -->  y

					m_min_i = std::min( i , m_min_i );
					m_min_x = std::min( x , m_min_x );
					m_min_y = std::min( y , m_min_y );
					m_max_i = std::max( i , m_max_i );
					m_max_x = std::max( x , m_max_x );
					m_max_y = std::max( y , m_max_y );
				}

	if(m_on_i.size() == 0) // grain does not exist
	{
		m_own_id = 0;
		return;
	}
	
	m_vis_i_top.resize(m_max_x - m_min_x + 1); BOOST_FOREACH(t_vis_slice& tmp,m_vis_i_top) tmp.resize(m_max_y - m_min_y + 1);// x,y  -->  i
	m_vis_x_top.resize(m_max_y - m_min_y + 1); BOOST_FOREACH(t_vis_slice& tmp,m_vis_x_top) tmp.resize(m_max_i - m_min_i + 1);// y,i  -->  x
	m_vis_y_top.resize(m_max_i - m_min_i + 1); BOOST_FOREACH(t_vis_slice& tmp,m_vis_y_top) tmp.resize(m_max_x - m_min_x + 1);// i,x  -->  y
	m_vis_i_bot.resize(m_max_x - m_min_x + 1); BOOST_FOREACH(t_vis_slice& tmp,m_vis_i_bot) tmp.resize(m_max_y - m_min_y + 1);// x,y  -->  i
	m_vis_x_bot.resize(m_max_y - m_min_y + 1); BOOST_FOREACH(t_vis_slice& tmp,m_vis_x_bot) tmp.resize(m_max_i - m_min_i + 1);// y,i  -->  x
	m_vis_y_bot.resize(m_max_i - m_min_i + 1); BOOST_FOREACH(t_vis_slice& tmp,m_vis_y_bot) tmp.resize(m_max_x - m_min_x + 1);// i,x  -->  y


//	for(t_coord i = m_min_i ; i <= m_max_i ; i++)
//		for(t_coord x = m_min_x ; x <= m_max_x ; x++)
//			for(t_coord y = m_min_y ; y <= m_max_y ; y++)

	for(t_coord i=0;i<si;i++) {
		for(t_coord x=0;x<sx;x++) {
			for(t_coord y=0;y<sy;y++) {
				if(dat[i][x][y] == m_own_id) {
					bool	i_top=false,i_bot=false,
						x_top=false,x_bot=false,
						y_top=false,y_bot=false;
					bool	i_topcon=false,i_botcon=false,
						x_topcon=false,x_botcon=false,
						y_topcon=false,y_botcon=false;

					if(i==0){// i
						i_bot=true;
					} else {
						if(dat[i-1][x][y] == 0) {
							i_bot=true;
						} else {
							if(dat[i-1][x][y] != dat[i][x][y])
								i_botcon=true;
//							if(dat[i-1][x][y] > dat[i][x][y])
//								note_contact_info(i,x,y,I_BOT,con);
						}
					}

					if(i==si-1){
						i_top=true;
					} else {
						if(dat[i+1][x][y] == 0) {
							i_top=true;
						} else {
							if(dat[i+1][x][y] != dat[i][x][y])
								i_topcon=true;
//							if(dat[i+1][x][y] > dat[i][x][y])
//								note_contact_info(i,x,y,I_TOP,con);
						}
					}// i

					if(x==0){// x
						x_bot=true;
					} else {
						if(dat[i][x-1][y] == 0) {
							x_bot=true;
						} else {
							if(dat[i][x-1][y] != dat[i][x][y])
								x_botcon=true;
//							if(dat[i][x-1][y] > dat[i][x][y])
//								note_contact_info(i,x,y,X_BOT,con);
						}
					}

					if(x==sx-1){
						x_top=true;
					} else {
						if(dat[i][x+1][y] == 0) {
							x_top=true;
						} else {
							if(dat[i][x+1][y] != dat[i][x][y])
								x_topcon=true;
//							if(dat[i][x+1][y] > dat[i][x][y])
//								note_contact_info(i,x,y,X_TOP,con);
						}
					}// x

					if(y==0){// y
						y_bot=true;
					} else {
						if(dat[i][x][y-1] == 0) {
							y_bot=true;
						} else {
							if(dat[i][x][y-1] != dat[i][x][y])
								y_botcon=true;
//							if(dat[i][x][y-1] > dat[i][x][y])
//								note_contact_info(i,x,y,Y_BOT,con);
						}
					}

					if(y==sy-1){
						y_top=true;
					} else {
						if(dat[i][x][y+1] == 0) {
							y_top=true;
						} else {
							if(dat[i][x][y+1] != dat[i][x][y])
								y_topcon=true;
//							if(dat[i][x][y+1] > dat[i][x][y])
//								note_contact_info(i,x,y,Y_TOP,con);
						}
					}// y

					if(i_top)    m_vis_i_top[ x - m_min_x ][ y - m_min_y ].push_back(vc(i, dat[i][x][y]   ));
					if(i_bot)    m_vis_i_bot[ x - m_min_x ][ y - m_min_y ].push_back(vc(i, dat[i][x][y]   ));
					if(x_top)    m_vis_x_top[ y - m_min_y ][ i - m_min_i ].push_back(vc(x, dat[i][x][y]   ));
					if(x_bot)    m_vis_x_bot[ y - m_min_y ][ i - m_min_i ].push_back(vc(x, dat[i][x][y]   ));
					if(y_top)    m_vis_y_top[ i - m_min_i ][ x - m_min_x ].push_back(vc(y, dat[i][x][y]   ));
					if(y_bot)    m_vis_y_bot[ i - m_min_i ][ x - m_min_x ].push_back(vc(y, dat[i][x][y]   ));

					if(i_topcon) m_vis_i_top[ x - m_min_x ][ y - m_min_y ].push_back(vc(i, dat[i+1][x][y] ));
					if(i_botcon) m_vis_i_bot[ x - m_min_x ][ y - m_min_y ].push_back(vc(i, dat[i-1][x][y] ));
					if(x_topcon) m_vis_x_top[ y - m_min_y ][ i - m_min_i ].push_back(vc(x, dat[i][x+1][y] ));
					if(x_botcon) m_vis_x_bot[ y - m_min_y ][ i - m_min_i ].push_back(vc(x, dat[i][x-1][y] ));
					if(y_topcon) m_vis_y_top[ i - m_min_i ][ x - m_min_x ].push_back(vc(y, dat[i][x][y+1] ));
					if(y_botcon) m_vis_y_bot[ i - m_min_i ][ x - m_min_x ].push_back(vc(y, dat[i][x][y-1] ));
				}
			}
		}
	}
	std::cerr << "grain surface extracted: " << ((int)(m_own_id)) << "\n";
}

void GrainSurface::draw_single_side(t_coord i,t_coord x,t_coord y,enum_side side,float d2)
{
	glTranslatef(x,y,i);
	glBegin(GL_QUADS);
	switch(side)
	{
		case I_TOP : glNormal3f(0,0, 1);glVertex3f(-d2,-d2, d2);glVertex3f( d2,-d2, d2);glVertex3f( d2, d2, d2);glVertex3f(-d2, d2, d2);break;
		case I_BOT : glNormal3f(0,0,-1);glVertex3f(-d2,-d2,-d2);glVertex3f( d2,-d2,-d2);glVertex3f( d2, d2,-d2);glVertex3f(-d2, d2,-d2);break;
		case X_TOP : glNormal3f( 1,0,0);glVertex3f( d2,-d2,-d2);glVertex3f( d2,-d2, d2);glVertex3f( d2, d2, d2);glVertex3f( d2, d2,-d2);break;
		case X_BOT : glNormal3f(-1,0,0);glVertex3f(-d2,-d2,-d2);glVertex3f(-d2,-d2, d2);glVertex3f(-d2, d2, d2);glVertex3f(-d2, d2,-d2);break;
		case Y_TOP : glNormal3f(0, 1,0);glVertex3f(-d2, d2,-d2);glVertex3f( d2, d2,-d2);glVertex3f( d2, d2, d2);glVertex3f(-d2, d2, d2);break;
		case Y_BOT : glNormal3f(0,-1,0);glVertex3f(-d2,-d2,-d2);glVertex3f( d2,-d2,-d2);glVertex3f( d2,-d2, d2);glVertex3f(-d2,-d2, d2);break;
	}
	glEnd();
	glTranslatef(-x,-y,-i);
}

void GrainSurface::plot_sides(t_visibility& vis,t_coord A,t_coord B,t_coord min_A,t_coord min_B,enum_side side,Config& config,SafeVectors3& colors)
{
	bool            names                 = config.draw_with_names();
	qglviewer::Vec  cam                   = config.view_direction();
	bool            plot_faster_with_back = config.draw_backfaces();
	short int       draw_what_surfaces    = config.draw_surface_or_contact_or_all();
	int             DD                    = config.display_steps();
	float           d2                    = DD*0.5;	
	int             CC(0);
	bool            bottom=false;
	switch(side)
	{
		case I_BOT : bottom=true;
		case I_TOP : CC=2;break;
		case X_BOT : bottom=true;
		case X_TOP : CC=0;break;
		case Y_BOT : bottom=true;
		case Y_TOP : CC=1;break;
	}
	for(t_coord a = 0 ; a<A ; a+=DD)
		for(t_coord b = 0 ; b < B ; b+=DD)
		{
			bool front = (bottom?(cam[CC]>0):(cam[CC]<0));
			if(plot_faster_with_back || front )
			{
				BOOST_FOREACH(vc& C, vis[a][b])
				{
					t_coord c     = C.level;
					t_voxel color = C.color;
					if (
							(draw_what_surfaces == 2) ||
							(draw_what_surfaces == 1 && color == m_own_id) ||
							(draw_what_surfaces == 0 && color != m_own_id)
					   )
					{
						if(!names)
						{
							Vector3r co=colors[color];//-Vector3r(1,1,1)*0.07;
							if(config.transparency())
								glColor4f(co[0],co[1],co[2],front?0.7:0.5);
							else
								glColor3f(co[0],co[1],co[2]);
							if(!config.display_selection())
								config.light_selection(m_own_id);
						}
						int i(0),x(0),y(0);
						switch(side)
						{
							case I_BOT : 
							case I_TOP : i=c;       x=a+min_A; y=b+min_B;break;
							case X_BOT : 
							case X_TOP : i=b+min_B; x=c;       y=a+min_A;break;
							case Y_BOT : 
							case Y_TOP : i=a+min_A; x=b+min_B; y=c;      break;
						}
						draw_single_side(i,x,y,side,d2);
					}
				}
			}
		}
}

void GrainSurface::draw(SafeVectors3& colors,Config& c)
{
	if(m_own_id == 0)
		return;

	plot_sides(m_vis_i_top,m_max_x - m_min_x + 1,m_max_y - m_min_y + 1,m_min_x,m_min_y,I_TOP,c,colors); // x,y  -->  i
	plot_sides(m_vis_i_bot,m_max_x - m_min_x + 1,m_max_y - m_min_y + 1,m_min_x,m_min_y,I_BOT,c,colors); // x,y  -->  i
	plot_sides(m_vis_x_top,m_max_y - m_min_y + 1,m_max_i - m_min_i + 1,m_min_y,m_min_i,X_TOP,c,colors); // y,i  -->  x
	plot_sides(m_vis_x_bot,m_max_y - m_min_y + 1,m_max_i - m_min_i + 1,m_min_y,m_min_i,X_BOT,c,colors); // y,i  -->  x
	plot_sides(m_vis_y_top,m_max_i - m_min_i + 1,m_max_x - m_min_x + 1,m_min_i,m_min_x,Y_TOP,c,colors); // i,x  -->  y
	plot_sides(m_vis_y_bot,m_max_i - m_min_i + 1,m_max_x - m_min_x + 1,m_min_i,m_min_x,Y_BOT,c,colors); // i,x  -->  y
}

//void note_contact_info(int i,int x,int y,E_SIDE side,T_CONTACTS& con)
//{
//	int other,current=chars[i][x][y];
//	Vector3r dir;
//	switch(side)
//	{
//		case I_BOT : other=chars[i-1][x][y];dir=Vector3r(0,0,-1);break;
//		case I_TOP : other=chars[i+1][x][y];dir=Vector3r(0,0, 1);break;
//		case X_BOT : other=chars[i][x-1][y];dir=Vector3r(-1,0,0);break;
//		case X_TOP : other=chars[i][x+1][y];dir=Vector3r( 1,0,0);break;
//		case Y_BOT : other=chars[i][x][y-1];dir=Vector3r(0,-1,0);break;
//		case Y_TOP : other=chars[i][x][y+1];dir=Vector3r(0, 1,0);break;
//	}
//	if(current>other)
//		std::cerr << "WARNING: note_contact_info: current>other, wtf?\n";
//
//	con[current][other].normal-=dir;
//	con[current][other].center+=Vector3r(x,y,i);
//	con[current][other].count+=1.0;
//	con[current][other].color=colors[other];//(colors[current]+colors[other])*0.5;
//
//
//	con[other][current].normal+=dir;
//	con[other][current].center+=Vector3r(x,y,i);
//	con[other][current].count+=1.0;
//	con[other][current].color=colors[current];
//}


YADE_REQUIRE_FEATURE(PHYSPAR);

