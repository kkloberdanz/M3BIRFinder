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

    // called from within startExecutables()
    void run_alignment(std::string);
    void run_full_align(std::string, std::string, std::string);
    void run_get_reads(std::string, std::string, std::string, std::string, bool);

  public:
    BWA();
    void startExecutables(int);
}; // class BWA
} // namespace bir
#endif //BWA_HPP
