#ifndef BWA_HPP
#define BWA_HPP

#include <mutex>

namespace bir {

class BWA {
  private:
    std::mutex mtx;
    std::vector<std::string> sReadsFile_v;

    int executeBwaIndex(std::string);
    int executeBwaAligner(std::string, std::string, std::string);
    int getReads(std::string, std::string, std::string, std::string, bool);
    int convertSAMtoFASTA(std::string);
    int filterOut (std::string, std::string, std::string, std::string);

  public:
    BWA();
    void startExecutables(int);
};
}
#endif //BWA_HPP
