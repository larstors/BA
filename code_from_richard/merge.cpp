// Merge two cluster files to stdout

#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char* argv[]) {

  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " cluster-file-1 cluster-file-2" << std::endl;
    return 1;
  }

  // Open the two files
  std::ifstream f1(argv[1]);
  if(!f1) {
    std::cerr << "Could not open '" << argv[1] << "' for reading" << std::endl;
    return 1;
  }
  std::ifstream f2(argv[2]);
  if(!f2) {
    std::cerr << "Could not open '" << argv[2] << "' for reading" << std::endl;
    return 1;
  }

  bool a1=true, a2=true; // Files are considered active after opening

  // Read a line at a time
  for(;;) {
    std::string l1, l2;
    a1 = bool(std::getline(f1, l1));
    if(!a1) break;
    a2 = bool(std::getline(f2, l2));
    if(!a2) {
      std::cout << l1 << std::endl;
      break;
    }

    // Check we have either (commment-comment) or (data-data)
    if((l1[0] == '#')!=(l2[0]=='#')) {
      std::cerr << "Location of comment lines differs" << std::endl;
      return 1;
    }
    if(l1[0] == '#') {
      // Comment - need to be identical
      if(l1!=l2) {
        std::cerr << "Content of comments differ" << std::endl;
        return 1;
      }
      // They are the same, echo
      std::cout << l1 << std::endl;
    } else {
      // Data
      std::istringstream s1(l1), s2(l2);
      bool sep = false;
      // Split the string into a sequence of integers; assume zero if one runs out before the other; output sum
      for(;;) {
        int i1, i2;
        s1 >> i1;
        if(s1.fail()) {
          std::cout << (s2.eof() ? "" : s2.str().substr(s2.tellg())) << std::endl;
          break;
        }
        s2 >> i2;
        if(s2.fail()) {
          std::cout << (sep ? " " : "") << i1;
          std::cout << (s1.eof() ? "" : s1.str().substr(s1.tellg())) << std::endl;
          break;
        }
        std::cout << (sep ? " " : "") << (i1+i2);
        sep = true;
      }
    }
  }

  // Write out any remaining lines unchanged
  while(a1||a2) {
    std::string l;
    if(a1 && (a1 = bool(std::getline(f1, l)))) std::cout << l << std::endl;
    if(a2 && (a2 = bool(std::getline(f2, l)))) std::cout << l << std::endl;
  }

  return 0;
}
