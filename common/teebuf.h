/**
 *  @file    teebuf.cc
 *  @author  Dietmar Kühl
 *  @date    2012-12-02
 *
 *  @brief Struct that can be used to write to streambuffers simultaneously.
 * 
 *  The expected usage would be this:
 * 
 *  // Create a teebuf from a ostringstream (or ofstream, etc.), then create a
 *  // ostream from this
 *  std::ostringstream sout;
 *  teebuf sbuf(sout.rdbuf(), std::cout.rdbuf());
 *  std::ostream out(&sbuf);
 *
 * // Save the original cout streambuf and assign the new ostream's streambuf
 * // to cout
 *  std::streambuf* old_cout_streambuf = std::cout.rdbuf();
 *  std::cout.rdbuf( out.rdbuf() );
 * 
 *  ... your code ...
 * 
 *  // Restore the original cout streambuf (otherwise you will segfault)
 *  std::cout.rdbuf(old_cout_streambuf);
 *
 */

#ifndef TEEBUF_H_
#define TEEBUF_H_

#include <streambuf>

struct teebuf : std::streambuf {
    std::streambuf* sb1_;
    std::streambuf* sb2_;

    teebuf(std::streambuf* sb1, std::streambuf* sb2) : sb1_(sb1), sb2_(sb2) {}
    int overflow(int c) {
        typedef std::streambuf::traits_type traits;
        bool rc(true);
        if (!traits::eq_int_type(traits::eof(), c)) {
            traits::eq_int_type(this->sb1_->sputc(c), traits::eof()) && (rc = false);
            traits::eq_int_type(this->sb2_->sputc(c), traits::eof()) && (rc = false);
        }
        return rc ? traits::not_eof(c) : traits::eof();
    }
    int sync() {
        bool rc(true);
        this->sb1_->pubsync() != -1 || (rc = false);
        this->sb2_->pubsync() != -1 || (rc = false);
        return rc ? 0 : -1;
    }
};

#endif  // TEEBUF_H_