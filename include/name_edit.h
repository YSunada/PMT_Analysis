#ifndef INCLUDED_NAME_EDIT_SUNA
#define INCLUDED_NAME_EDIT_SUNA

#include "TString.h"
#include "string"

///////////////////////////////////
//     Get Serial Number         //
///////////////////////////////////
TString GetSerial1(TString file_name)
{
  std::string str_filename = file_name.Data();
  int num = str_filename.rfind(".root");
  TString Serial = str_filename.substr(num - 16,6);
  return Serial;
}

TString GetSerial1(std::string file_name)
{
  int num = file_name.rfind(".root");
  TString Serial = file_name.substr(num - 16,6);
   return Serial;
}

TString GetSerial2(TString file_name)
{
  std::string str_filename = file_name.Data();
  int num = str_filename.rfind(".root");
  TString Serial = str_filename.substr(num - 7,6);
  return Serial;
}

TString GetSerial2(std::string file_name)
{
  int num = file_name.rfind(".root");
  TString Serial = file_name.substr(num - 7,6);
  return Serial;
}

TString GetSerial(std::string file_name)
{
  int num = file_name.rfind(".root");
  TString Serial = file_name.substr(num - 18,6);
   return Serial;
}

TString GetSerial(TString file_name)
{
  std::string str_filename = file_name.Data();
  int num = str_filename.rfind(".root");
  TString Serial = str_filename.substr(num - 19,6);
  return Serial;
}


///////////////////////////////////
//         Get Date              //
///////////////////////////////////

TString GetDate(TString file_name)
{
  std::string str_filename = file_name.Data();
  int num = str_filename.rfind(".root");
  TString Serial = str_filename.substr(num - 28,12);
  return Serial;
}

TString GetDate(std::string file_name)
{
  int num = file_name.rfind(".root");
  TString Serial = file_name.substr(num - 28,12);
  return Serial;
}

TString GetDate2(TString file_name)
{
  std::string str_filename = file_name.Data();
  int num = str_filename.rfind(".root");
  TString Serial = str_filename.substr(num - 12,12);
  return Serial;
}

///////////////////////////////////
//       Get Data Type           //
///////////////////////////////////
TString GetDataType(TString file_name)
{
  std::string str_filename = file_name.Data();
  int num = str_filename.rfind(".root");
  if(file_name[num-1] == 'g')
    {
      return "OG";
    }
  else if(file_name[num-1] == 'a')
    {
      return "AP";
    }
  else
    {
      return "Unknown";
    }
}


TString GetDataType(std::string file_name)
{
  int num = file_name.rfind(".root");
  if(file_name[num-1] == 'g')
    {
      return "OG";
    }
  else if(file_name[num-1] == 'a')
    {
      return "AP";
    }
  else
    {
      return file_name[num-1];
    }
}

void GetSerialAndDate(TString *serial,TString *date,int ch,TString filename)
{
if(ch == 0)
	    {
	      *serial = GetSerial1(filename);
	    }
	  else if(ch == 1)
	    {
	      *serial = GetSerial2(filename);
	    }
	  *date = GetDate(filename);  
}

#endif
