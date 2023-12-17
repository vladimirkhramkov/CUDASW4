#ifndef READER_INTERFACE_H
#define READER_INTERFACE_H

#include <string>

class ReaderInterface {
public:
    virtual ~ReaderInterface() {}
    virtual int next() = 0;
    virtual std::string& getCurrentHeader() = 0;
    virtual std::string& getCurrentSequence() = 0;
};

#endif // READER_INTERFACE_H
