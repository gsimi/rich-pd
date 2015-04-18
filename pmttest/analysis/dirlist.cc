
/* helper functions to list the contents of a directory */

Bool_t IsItDirectory(const char *name, char * dirfile) 
{
  // Check if name is a directory.
  Long64_t size;
  Long_t id, flags, modtime;

  gSystem->ChangeDirectory(dirfile);
  flags = id = size = modtime = 0;
  gSystem->GetPathInfo(name, &id, &size, &flags, &modtime);
  Int_t isdir = (Int_t)flags & 2;

  return isdir ? kTRUE : kFALSE;
}
 
vector<string> GetListOfFiles(char* dirname, char* match=0){
  string pwd(gSystem->pwd());
  void *dir = gSystem->OpenDirectory(dirname);
  if (!dir) {
    vector<string> empty;
    return  empty;
  }

  const char *file = 0;
  vector<string> contents;
  while ((file = gSystem->GetDirEntry(dir))) {
    if (!( IsItDirectory(file,dirname)) ) {
      string filestr(file);
      if ((match==0) || (filestr.find(match)!=string::npos) ){
	string filepath(dirname);
	filepath.append(file);
	contents.push_back(filepath);
      }
    }
  }
  gSystem->FreeDirectory(dir);
  gSystem->ChangeDirectory(pwd.c_str());
  return contents;
}

