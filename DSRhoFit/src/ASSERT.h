#ifndef ASSERT_H_INCLUDED
    #define ASSERT_H_INCLUDED
    #ifndef DEBUG
        #define ASSERT(x)
    #else
        #include <iostream>
        #define ASSERT(x) \
        if (! (x)) \
        { \
            std::cout << "ERROR!! Assert " << #x << " failed\n"; \
            std::cout << " on line " << __LINE__  << "\n"; \
            std::cout << " in file " << __FILE__ << "\n";  \
        }
    #endif
#endif
