// FastaLine.cc                                                   -*-C++-*-

#include "linereader.hh"
#include "utilities.hh"

#include <cassert>
#include <iostream>
#include <iomanip>


std::ostream& Line::print(std::ostream& stream) const
{
    if (stream.good()) {
        stream << this->line << std::endl;
       // stream << orfid << "\t" << line << "\n";
       // stream << line << "\n";
    }

	return stream;
}

// FREE OPERATORS
std::istream& operator>>(std::istream& stream, Line& rhs)
{
    char buf[10000];
    std::string _line;
    if (getline(stream, _line)) {
        //rhs.line.assign(_line.begin(), _line.end());
        rhs.line = _line;
        char *field = split_n_pick(_line, buf, '\t', 0);
        rhs.orfid = string(field);
        //std::cout << _line << std::endl;
    }
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Line& rhs)
{
    stream << rhs.line;
    return stream;
}

 // close package namespace
