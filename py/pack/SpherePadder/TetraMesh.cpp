/*************************************************************************
 *  Copyright (C) 2009 by Jean Francois Jerier                            *
 *  jerier@geo.hmg.inpg.fr                                                *
 *  Copyright (C) 2009 by Vincent Richefeu                                *
 *  vincent.richefeu@hmg.inpg.fr                                          *
 *                                                                        *
 *  This program is free software; it is licensed under the terms of the  *
 *  GNU General Public License v2 or later. See file LICENSE for details. *
 *************************************************************************/

#include "TetraMesh.hpp"

TetraMesh::TetraMesh ()
{
    isOrganized = false;
}

// FIXME Some problem should occur with gmsh because nodes are sometimes duplicated.
void TetraMesh::read_gmsh (const char* name)
{
    ifstream meshFile(name);
    if(!meshFile)
    {
        cerr << "TetraMesh::read_gmsh, cannot open file " << name << endl;
        return;
    }

    string token;
    char not_read[150];
    meshFile >> token;

    while(meshFile)
    {
        if (token == "$Nodes")
        {
            unsigned int nbnodes;
            unsigned int num_node;
            Node N;

            meshFile >> nbnodes;
            for (unsigned int n = 0 ; n < nbnodes ; ++n)
            {
                meshFile >> num_node >> N.x >> N.y >> N.z;
                // cout << num_node << " " << N.x << "" << N.y << " " << N.z << endl;
                node.push_back(N);
            }
        }

        if (token == "$Elements")
        {
            unsigned int nbElements;
            unsigned int num_element, element_type, nbTags ;
            Tetraedre T;
            unsigned int t = 0;

            meshFile >> nbElements;
            for (unsigned int e = 0 ; e < nbElements ; ++e)
            {
                meshFile >> num_element >> element_type;
                // 4-node tetrahedron
                if (element_type != 4)
                {
                    meshFile.getline(not_read,150);
                    continue;
                }

                meshFile >> nbTags;
                // the third tag is the number of a mesh partition to which the element belongs
                unsigned int tag;
                for (unsigned int tg = 0 ; tg < nbTags ; ++(tg))
                    { meshFile >> tag; }

                    meshFile >> T.nodeId[0] >> T.nodeId[1] >> T.nodeId[2] >> T.nodeId[3];

                // nodeId has 0-offset
                // (0 in C/C++ corresponds to 1 in the file)
                T.nodeId[0] -= 1;
                T.nodeId[1] -= 1;
                T.nodeId[2] -= 1;
                T.nodeId[3] -= 1;

                node[T.nodeId[0]].tetraOwner.push_back(t);
                node[T.nodeId[1]].tetraOwner.push_back(t);
                node[T.nodeId[2]].tetraOwner.push_back(t);
                node[T.nodeId[3]].tetraOwner.push_back(t);

                tetraedre.push_back(T);
                ++t;
            }
        }

        if (token == "$EndElements") break;

        meshFile >> token;
    }

    organize ();
}


void TetraMesh::read_inp (const char* name)
{
  ifstream meshFile(name);
  if(!meshFile)
  {
	cerr << "TetraMesh::read_inp, cannot open file " << name << endl;
	return;
  }
  
  char line[256];
  
  meshFile.getline(line,256);
  while (!meshFile.eof() && line[0] == '*') meshFile.getline(line,256);
  
  Node N;
  unsigned int ID;
  while (!meshFile.eof() && line[0] != '*')
  {
	sscanf(line,"%d, %lf, %lf, %lf",&ID,&N.x,&N.y,&N.z);
	node.push_back(N);
	meshFile.getline(line,256);
  }

  while (!meshFile.eof() && line[0] == '*') meshFile.getline(line,256);

  Tetraedre T;
  while (!meshFile.eof() && line[0] != '*')
  {
	sscanf(line,"%d, %d, %d, %d, %d",&ID,&T.nodeId[0],&T.nodeId[1],&T.nodeId[2],&T.nodeId[3]);

	// nodeId has 0-offset
	// (0 in C/C++ corresponds to 1 in the file)
	T.nodeId[0] -= 1;
	T.nodeId[1] -= 1;
	T.nodeId[2] -= 1;
	T.nodeId[3] -= 1;

	ID = ID-1;
	node[T.nodeId[0]].tetraOwner.push_back(ID);
	node[T.nodeId[1]].tetraOwner.push_back(ID);
	node[T.nodeId[2]].tetraOwner.push_back(ID);
	node[T.nodeId[3]].tetraOwner.push_back(ID);

	tetraedre.push_back(T);
	meshFile.getline(line,256);
  }
  
  organize ();
}


void TetraMesh::read (const char* name)
{
    ifstream meshFile(name);
    if(!meshFile)
    {
        cerr << "TetraMesh::read_data, cannot open file " << name << endl;
        return;
    }

    string token;
    meshFile >> token;
    bool with_numbers = false;
    unsigned int number;

    while(meshFile)
    {
        if (token == "WITH_NUMBERS") with_numbers = true;
        if (token == "NODES")
        {
            unsigned int nbnodes;
            Node N;

            meshFile >> nbnodes;
            for (unsigned int n = 0 ; n < nbnodes ; ++n)
            {
                if (with_numbers) meshFile >> number;
                meshFile >> N.x >> N.y >> N.z;

                //cout << n << " " << N.x << " " << N.y << " " << N.z << endl;
                node.push_back(N);
            }
        }

        if (token == "TETRA")
        {
            unsigned int nbTetra;
            unsigned int nbignored = 0;
            Tetraedre T;

            meshFile >> nbTetra;
            for (unsigned int t = 0 ; t < nbTetra ; ++t)
            {
                if (with_numbers) meshFile >> number;
                meshFile >> T.nodeId[0] >> T.nodeId[1] >> T.nodeId[2] >> T.nodeId[3];

                //cout << t << " | " << T.nodeId[0] << " " <<  T.nodeId[1] << " " << T.nodeId[2] << " " << T.nodeId[3] << endl;

                // numbers begin at 0 instead of 1
                // (0 in C/C++ corresponds to 1 in the file)
                T.nodeId[0] -= 1;
                T.nodeId[1] -= 1;
                T.nodeId[2] -= 1;
                T.nodeId[3] -= 1;

                if (T.nodeId[0] >= node.size() || T.nodeId[1] >= node.size() || T.nodeId[2] >= node.size() || T.nodeId[3] >= node.size())
                {
                    nbignored++;
                    //cerr << "(1) error in mesh file, tetrahedron ignored! nbignored = " << nbignored << endl;
                    continue;
                }
                if (T.nodeId[0] < 0 || T.nodeId[1] < 0 || T.nodeId[2] < 0 || T.nodeId[3] < 0)
                {
                    nbignored++;
                    //cerr << "(2) error in mesh file, tetrahedron ignored! nbignored = " << nbignored << endl;
                    continue;
                }
                //cout << "num = " << (t-nbignored) << endl;

                node[T.nodeId[0]].tetraOwner.push_back(t-nbignored);
                node[T.nodeId[1]].tetraOwner.push_back(t-nbignored);
                node[T.nodeId[2]].tetraOwner.push_back(t-nbignored);
                node[T.nodeId[3]].tetraOwner.push_back(t-nbignored);

                tetraedre.push_back(T);

            }
        }

        if (token == "EOF") break;

        meshFile >> token;
    }

    organize ();
}


void TetraMesh::write_surface_MGP (const char* name)
{
    ofstream fmgpost(name);

    fmgpost << "<?xml version=\"1.0\"?>" << endl
        << " <mgpost mode=\"3D\">" << endl
        << "  <state id=\"" << 1
        << "\" time=\"" << 0.0 << "\">" << endl;

    double xm,ym,zm;
    unsigned int id1,id2,id3;
    unsigned int idPolyhedron = 1;
    double div3 = 0.3333333333333333;
    for (unsigned int i = 0 ; i < face.size() ; ++i)
    {
        if(!(face[i].belongBoundary)) continue;
        id1 = face[i].nodeId[0];
        id2 = face[i].nodeId[1];
        id3 = face[i].nodeId[2];
        xm = div3 * (node[id1].x + node[id2].x + node[id3].x);
        ym = div3 * (node[id1].y + node[id2].y + node[id3].y);
        zm = div3 * (node[id1].z + node[id2].z + node[id3].z);

        fmgpost << "   <body>" << endl;
        fmgpost << "    <POLYE id=\"" << idPolyhedron++ << "\" r=\"" << min_segment_length << "\">" << endl
            << "     <position x=\"" << xm + xtrans << "\" y=\"" << ym + ytrans << "\" z=\"" << zm + ztrans << "\"/>" << endl
            << "     <node x=\"" << node[id1].x - xm << "\" y=\"" << node[id1].y - ym << "\" z=\"" << node[id1].z - zm << "\"/>" << endl
            << "     <node x=\"" << node[id2].x - xm << "\" y=\"" << node[id2].y - ym << "\" z=\"" << node[id2].z - zm << "\"/>" << endl
            << "     <node x=\"" << node[id3].x - xm << "\" y=\"" << node[id3].y - ym << "\" z=\"" << node[id3].z - zm << "\"/>" << endl;

        if (face[i].normal_swap)
            fmgpost << "     <face n1=\"" << 0 << "\" n2=\"" << 2 << "\" n3=\"" << 1 << "\"/>" << endl;
        else
            fmgpost << "     <face n1=\"" << 0 << "\" n2=\"" << 1 << "\" n3=\"" << 2 << "\"/>" << endl;

        fmgpost << "    </POLYE>" << endl << flush;

        fmgpost << "   </body>" << endl;
    }

    fmgpost << "  </state>" << endl
        << " </mgpost>" << endl;
}


int compareInt (const void * a, const void * b)
{
    return ( *(int*)a > *(int*)b ) ? 1 :-1;
}


void TetraMesh::organize ()
{
	// Display information
	cout << "nb Nodes = " << node.size() << endl;
	cout << "nb Tetra = " << tetraedre.size() << endl;

	if(node.size()==0) throw std::runtime_error("TetraMesh has 0 vertices (undetected error importing mesh file?)");
  
    cout << "Translate all nodes in such a manner that all coordinates are > 0" << endl;
	// (Don't know if it is absolutly necessary)
    xtrans = node[0].x;
    ytrans = node[0].y;
    ztrans = node[0].z;
    for (unsigned int i = 1 ; i < node.size() ; ++i)
    {
        xtrans = (xtrans < node[i].x) ? xtrans : node[i].x;
        ytrans = (ytrans < node[i].y) ? ytrans : node[i].y;
        ztrans = (ztrans < node[i].z) ? ztrans : node[i].z;
    }
    xtrans = (xtrans < 0.0) ? -xtrans : 0.0;
    ytrans = (ytrans < 0.0) ? -ytrans : 0.0;
    ztrans = (ztrans < 0.0) ? -ztrans : 0.0;
    for (unsigned int i = 0 ; i < node.size() ; ++i)
    {
        node[i].x += xtrans;
        node[i].y += ytrans;
        node[i].z += ztrans;
    }

    cout << "Organize tetraedre nodes in ascending order" << endl;
    for (unsigned int i = 0 ; i < tetraedre.size() ; ++i)
    {
        qsort (&(tetraedre[i].nodeId[0]), 4, sizeof(int), compareInt);
    }

    // Face creation
    vector <Face> tmpFace;       // This will contain all faces more than one time (with duplications)
    Face F;
    F.tetraOwner.push_back(0);
    F.belongBoundary = true;
    F.normal_swap    = false;
    for (unsigned int i = 0 ; i < tetraedre.size() ; ++i)
    {
        F.tetraOwner[0] = i;

        // Notice that nodes are still organized in ascending order
        F.nodeId[0] = tetraedre[i].nodeId[0];
        F.nodeId[1] = tetraedre[i].nodeId[1];
        F.nodeId[2] = tetraedre[i].nodeId[2];
        tmpFace.push_back(F);

        F.nodeId[0] = tetraedre[i].nodeId[1];
        F.nodeId[1] = tetraedre[i].nodeId[2];
        F.nodeId[2] = tetraedre[i].nodeId[3];
        tmpFace.push_back(F);

        F.nodeId[0] = tetraedre[i].nodeId[0];
        F.nodeId[1] = tetraedre[i].nodeId[2];
        F.nodeId[2] = tetraedre[i].nodeId[3];
        tmpFace.push_back(F);

        F.nodeId[0] = tetraedre[i].nodeId[0];
        F.nodeId[1] = tetraedre[i].nodeId[1];
        F.nodeId[2] = tetraedre[i].nodeId[3];
        tmpFace.push_back(F);
    }


// FIXME the following double loop is very slow for large meshes (more than 200 000 tetrahedrons)
// utiliser le set<Face,cmp_face>
	cout << "Search for duplicated and boundary faces" << endl;
    face.push_back(tmpFace[0]);
    bool duplicated;
    for (unsigned int i = 1 ; i < tmpFace.size() ; ++i)
    {
        duplicated = false;
        for (unsigned int n = 0 ; n < face.size() ; ++n)
        {
            if (   tmpFace[i].nodeId[0] == face[n].nodeId[0]
                && tmpFace[i].nodeId[1] == face[n].nodeId[1]
                && tmpFace[i].nodeId[2] == face[n].nodeId[2])
            {
                duplicated = true;
                face[n].tetraOwner.push_back(tmpFace[i].tetraOwner[0]);
                face[n].belongBoundary = false;
                break;
            }
        }

        if (!duplicated)
        {
            face.push_back(tmpFace[i]);
        }
    }
    tmpFace.clear();  // It should be useful to save memory before the end of the function

    cout <<  "Computation of unit vectors normal to the boundary faces" << endl;
    unsigned int n1,n2,n3,n4;
    unsigned int s1,s2,s3,s4;
    Tetraedre current_tetra;
    double v12[3], v13[3], v14[3], vect_product[3];
    double scal_product;
    for (unsigned int f = 0 ; f < face.size() ; ++f)
    {
        if (!(face[f].belongBoundary)) continue;

        current_tetra = tetraedre[ face[f].tetraOwner[0] ];

        n1 = current_tetra.nodeId[0];
        n2 = current_tetra.nodeId[1];
        n3 = current_tetra.nodeId[2];
        n4 = current_tetra.nodeId[3];

        s1 = face[f].nodeId[0];
        s2 = face[f].nodeId[1];
        s3 = face[f].nodeId[2];

        if      ((n1 != s1) && (n1 != s2) && (n1 != s3)) {s4 = n1;}
        else if ((n2 != s1) && (n2 != s2) && (n2 != s3)) {s4 = n2;}
        else if ((n3 != s1) && (n3 != s2) && (n3 != s3)) {s4 = n3;}
        else                                             {s4 = n4;}

        v12[0] = node[s2].x - node[s1].x;
        v12[1] = node[s2].y - node[s1].y;
        v12[2] = node[s2].z - node[s1].z;

        v13[0] = node[s3].x - node[s1].x;
        v13[1] = node[s3].y - node[s1].y;
        v13[2] = node[s3].z - node[s1].z;

        v14[0] = node[s4].x - node[s1].x;
        v14[1] = node[s4].y - node[s1].y;
        v14[2] = node[s4].z - node[s1].z;

        vect_product[0] = v12[1]*v13[2] - v12[2]*v13[1];
        vect_product[1] = v12[2]*v13[0] - v12[0]*v13[2];
        vect_product[2] = v12[0]*v13[1] - v12[1]*v13[0];

        scal_product = v14[0]*vect_product[0] + v14[1]*vect_product[1] + v14[2]*vect_product[2];

        if (scal_product > 0.0) face[f].normal_swap = true;
    }

    // Define for each node the faces it belong to
    for (unsigned int f = 0 ; f < face.size() ; ++f)
    {
        node[face[f].nodeId[0]].faceOwner.push_back(f);
        node[face[f].nodeId[1]].faceOwner.push_back(f);
        node[face[f].nodeId[2]].faceOwner.push_back(f);
    }

    // Segment creation
    vector <Segment> tmpSegment;
    Segment S;
    S.faceOwner.push_back(0);

    for (unsigned int i = 0 ; i < face.size() ; ++i)
    {
        S.faceOwner[0] = i;
        S.nodeId[0] = face[i].nodeId[0];
        S.nodeId[1] = face[i].nodeId[1];
        tmpSegment.push_back(S);

        S.nodeId[0] = face[i].nodeId[1];
        S.nodeId[1] = face[i].nodeId[2];
        tmpSegment.push_back(S);

        S.nodeId[0] = face[i].nodeId[0];
        S.nodeId[1] = face[i].nodeId[2];
        tmpSegment.push_back(S);
    }

// FIXME double loop can be very slow !!!!!
	cout << "Search for duplicated segments" << endl;
    segment.push_back(tmpSegment[0]);
    for (unsigned int i = 1 ; i < tmpSegment.size() ; ++i)
    {
        duplicated = false;
        for (unsigned int n = 0 ; n < segment.size() ; ++n)
        {
            if (   tmpSegment[i].nodeId[0] == segment[n].nodeId[0]
                && tmpSegment[i].nodeId[1] == segment[n].nodeId[1])
            {
                duplicated = true;
                segment[n].faceOwner.push_back(tmpSegment[i].faceOwner[0]);
                break;
            }

        }

        if (!duplicated)
        {
            segment.push_back(tmpSegment[i]);
        }
    }

    for (unsigned int s = 0 ; s < segment.size() ; ++s)
    {
        node[segment[s].nodeId[0]].segmentOwner.push_back(s);
        node[segment[s].nodeId[1]].segmentOwner.push_back(s);
    }

    cout << "Compute length of segments" << endl;
    double lx,ly,lz;
    unsigned int id1,id2;

    // length_min = length(0)
    id1 = segment[0].nodeId[0];
    id2 = segment[0].nodeId[1];
    lx  = node[id2].x - node[id1].x;
    ly  = node[id2].y - node[id1].y;
    lz  = node[id2].z - node[id1].z;
    min_segment_length = sqrt(lx*lx + ly*ly + lz*lz);

    mean_segment_length = 0.0;
    max_segment_length = 0.0;

    for (unsigned int s = 0 ; s < segment.size() ; ++s)
    {
        id1 = segment[s].nodeId[0];
        id2 = segment[s].nodeId[1];

        lx  = node[id2].x - node[id1].x;
        ly  = node[id2].y - node[id1].y;
        lz  = node[id2].z - node[id1].z;

        segment[s].length = sqrt(lx*lx + ly*ly + lz*lz);

        mean_segment_length += segment[s].length;
        min_segment_length = (segment[s].length < min_segment_length) ? segment[s].length : min_segment_length;
        max_segment_length = (segment[s].length > max_segment_length) ? segment[s].length : max_segment_length;
    }
    mean_segment_length /= (double)(segment.size());

//     cerr << "mean_segment_length = " << mean_segment_length << endl;
//     cerr << "min_segment_length = " << min_segment_length << endl;
//     cerr << "max_segment_length = " << max_segment_length << endl;

    // Define tetraedre neighbors FIXME still useful??
    bool stop = false;
    for (unsigned int t1 = 0 ; t1 < tetraedre.size() ; ++t1)
    {
        tetraedre[t1].tetraNeighbor.push_back(t1);
        for (unsigned int t2 = 0 ; t2 < tetraedre.size() ; ++t2)
        {
            /*
            if (   (tetraedre[t1].nodeId[0] > tetraedre[t2].nodeId[3])
            || (tetraedre[t1].nodeId[3] < tetraedre[t2].nodeId[0]) ) continue;
                  */

            // FIXME mettre du while... (?)
            stop = false;
            for (unsigned int i = 0 ; i < 4 ; i++)
            {
                for (unsigned int j = 0 ; j < 4 ; j++)
                {
                    if (tetraedre[t1].nodeId[i] == tetraedre[t2].nodeId[j])
                    {
                        if (t1 != t2) tetraedre[t1].tetraNeighbor.push_back(t2);
                        //if (t1 != t2) tetraedre[t2].tetraNeighbor.push_back(t1);
                        stop = true;
                        break;
                    }
                }
                if (stop) break;
            }
        }
    }

    /*
            for (unsigned int t1 = 0 ; t1 < tetraedre.size() ; ++t1)
            {
              cerr << "Tetra " << t1 << " a " <<  tetraedre[t1].tetraNeighbor.size() << " voisins" << endl;
              for (unsigned int i = 0 ; i < tetraedre[t1].tetraNeighbor.size() ; i++)
              {
                cerr <<  tetraedre[t1].tetraNeighbor[i] << ' ';
              }
              cerr << endl;
            }
            exit(0);
    */

    isOrganized = true;
}
