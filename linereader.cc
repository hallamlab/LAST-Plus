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
    }

	return stream;
}

// FREE OPERATORS
std::istream& operator>>(std::istream& stream, Line& rhs)
{
    char buf[10000];
    std::string _line;
    if (getline(stream, _line)) {
        rhs.line = _line;
        char *field = split_n_pick(_line, buf, '\t', 0);
        rhs.orfid = string(field);
        rhs.evalue = evalue_extractor_from_blast(_line);
    }
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Line& rhs)
{
    stream << rhs.line;
    return stream;
}

 // close package namespace
