// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#include "Wm3FoundationPCH.h"
#include "Wm3String.h"
using namespace Wm3;

//----------------------------------------------------------------------------
String::String (const char* acText)
{
    if (acText)
    {
        m_iLength = (int)strlen(acText);
        m_acText = WM3_NEW char[m_iLength+1];
        const size_t uiSize = (size_t)(m_iLength + 1);
        System::Strcpy(m_acText,uiSize,acText);
    }
    else
    {
        m_iLength = 0;
        m_acText = WM3_NEW char[1];
        m_acText[0] = 0;
    }
}
//----------------------------------------------------------------------------
String::String (int iLength, const char* acText)
{
    if (iLength < 0 || (iLength > 0 && !acText))
    {
        assert(false);
        iLength = 0;
    }

#ifdef _DEBUG
    // make sure acText has as many non-null characters as specified
    if (iLength > 0)
    {
        int i;
        for (i = 0; i < iLength; i++)
        {
            if (acText[i] == 0)
            {
                break;
            }
        }
        assert(i == iLength);
    }
#endif

    m_iLength = iLength;
    if (m_iLength > 0)
    {
        const size_t uiSrcSize = (size_t)m_iLength;
        const size_t uiDstSize = (size_t)(m_iLength + 1);
        m_acText = WM3_NEW char[uiDstSize];
        System::Strncpy(m_acText,uiDstSize,acText,uiSrcSize);
        m_acText[m_iLength] = 0;
    }
    else
    {
        m_acText = WM3_NEW char[1];
        m_acText[0] = 0;
    }
}
//----------------------------------------------------------------------------
String::String (const String& rkString)
{
    m_acText = 0;
    *this = rkString;
}
//----------------------------------------------------------------------------
String::~String ()
{
    WM3_DELETE[] m_acText;
}
//----------------------------------------------------------------------------
String& String::operator= (const String& rkString)
{
    m_iLength = rkString.m_iLength;
    WM3_DELETE[] m_acText;
    const size_t uiSize = (size_t)(m_iLength + 1);
    m_acText = WM3_NEW char[uiSize];
    System::Strcpy(m_acText,uiSize,rkString.m_acText);
    return *this;
}
//----------------------------------------------------------------------------
String& String::operator+= (const String& rkString)
{
    m_iLength += rkString.m_iLength;
    const size_t uiSize = (size_t)(m_iLength + 1);
    char* acNew = WM3_NEW char[uiSize];
    System::Strcpy(acNew,uiSize,m_acText);
    WM3_DELETE[] m_acText;
    m_acText = System::Strcat(acNew,uiSize,rkString.m_acText);   
    return *this;
}
//----------------------------------------------------------------------------
String String::operator+ (const String& rkString)
{
    int iLength = m_iLength + rkString.m_iLength;
    const size_t uiSize = (size_t)(iLength + 1);
    char* acNew = WM3_NEW char[uiSize];
    System::Strcpy(acNew,uiSize,m_acText);
    System::Strcat(acNew,uiSize,rkString.m_acText);
    String kAppended(iLength,acNew);
    WM3_DELETE[] acNew;
    return kAppended;
}
//----------------------------------------------------------------------------
bool String::operator== (const String& rkString) const
{
    return m_iLength == rkString.m_iLength
        && memcmp(m_acText,rkString.m_acText,m_iLength) == 0;
}
//----------------------------------------------------------------------------
bool String::operator!= (const String& rkString) const
{
    return m_iLength != rkString.m_iLength
        || memcmp(m_acText,rkString.m_acText,m_iLength) != 0;
}
//----------------------------------------------------------------------------
String::operator unsigned int () const
{
    unsigned int uiKey = 0;
    for (int i = 0; i < m_iLength; i++)
    {
        uiKey += (uiKey << 5) + m_acText[i];
    }
    return uiKey;
}
//----------------------------------------------------------------------------
String String::ToUpper () const
{
    String kCopy(*this);
    for (int i = 0; i < m_iLength; i++)
    {
        kCopy.m_acText[i] = toupper(m_acText[i]);
    }
    return kCopy;
}
//----------------------------------------------------------------------------
String String::ToLower () const
{
    String kCopy(*this);
    for (int i = 0; i < m_iLength; i++)
    {
        kCopy.m_acText[i] = tolower(m_acText[i]);
    }
    return kCopy;
}
//----------------------------------------------------------------------------

