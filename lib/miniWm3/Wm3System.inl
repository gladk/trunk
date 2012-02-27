// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

//----------------------------------------------------------------------------
template <class T>
void Allocate (int iCols, int iRows, T**& raatArray)
{
    raatArray = WM3_NEW T*[iRows];
    raatArray[0] = WM3_NEW T[iRows*iCols];
    for (int iRow = 1; iRow < iRows; iRow++)
    {
        raatArray[iRow] = &raatArray[0][iCols*iRow];
    }
}
//----------------------------------------------------------------------------
template <class T>
void Deallocate (T**& raatArray)
{
    WM3_DELETE[] raatArray[0];
    WM3_DELETE[] raatArray;
    raatArray = 0;
}
//----------------------------------------------------------------------------
template <class T>
void Allocate (int iCols, int iRows, int iSlices, T***& raaatArray)
{
    raaatArray = WM3_NEW T**[iSlices];
    raaatArray[0] = WM3_NEW T*[iSlices*iRows];
    raaatArray[0][0] = WM3_NEW T[iSlices*iRows*iCols];
    for (int iSlice = 0; iSlice < iSlices; iSlice++)
    {
        raaatArray[iSlice] = &raaatArray[0][iRows*iSlice];
        for (int iRow = 0; iRow < iRows; iRow++)
        {
            raaatArray[iSlice][iRow] =
                &raaatArray[0][0][iCols*(iRow+iRows*iSlice)];
        }
    }
}
//----------------------------------------------------------------------------
template <class T>
void Deallocate (T***& raaatArray)
{
    WM3_DELETE[] raaatArray[0][0];
    WM3_DELETE[] raaatArray[0];
    WM3_DELETE[] raaatArray;
    raaatArray = 0;
}
//----------------------------------------------------------------------------

