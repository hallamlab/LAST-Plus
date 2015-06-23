// FastaLine.hh                                                    -*-C++-*-

#ifndef INCLUDED_FASTA_Line
#define INCLUDED_FASTA_Line
#include <string>
#include <iterator>
#include "utilities.hh"

using namespace std;

class Line {
    // This attribute class...

    // DATA

    // FRIENDS
    friend std::istream& operator>>(std::istream&, Line&);

  public:
    std::string orfid;      // name of the Line
    std::string line;  // sequence of 'A', 'C', 'G', and 'T' characters
		double evalue;
    // CREATORS
    explicit Line();
        // Create a new 'Line' object having the (default) attribute values:
        //..
        //  name()     == "";
        //  sequence() == "";
        //..

    Line(const std::string& orfid, const std::string& line);
        // Create a new 'Line' object having the specified 'name' and
        // 'sequence' attribute values.

    Line(const Line& original);
        // Create a new 'Line' object having the same value as the specified
        // 'original' object.

    //! ~Line() = default;
        // Destroy this object.

    // MANIPULATORS
    Line& operator=(const Line& rhs);
        // Set the value of this object to the value of the specified 'rhs'
        // object, and return a reference providing modifiable access to this
        // object.

    void setOrfId(const std::string& value);
        // Set the 'name' attribute of this object to the specified 'value'.

    void setLine(const std::string& value);
        // Set the 'sequence' attribute of this object to the specified
        // 'value'.

		void setEvalue(const double& value);

    // ACCESSORS
    const std::string& getOrfId() const;
        // Return a reference providing non-modifiable access to the 'name'
        // attribute of this object.

    const std::string& getLine() const;
        // Return a reference providing non-modifiable access to the 'sequence'
        // attribute of this object.

                                  // Aspects

    /*std::ostream& print(std::ostream& stream,
                        int           level = 0,
                        int           spacesPerLevel = 4) const;*/

    std::ostream& print(std::ostream& stream) const;

};

// FREE OPERATORS
bool operator==(const Line& lhs, const Line& rhs);
    // Return 'true' if the specified 'lhs' and 'rhs' objects have the same
    // value and 'false' otherwise.  Two 'Line' objects have the same value
    // if the corresponding values of their 'sequence' and 'name' attributes have
    // the same value.

bool operator!=(const Line& lhs, const Line& rhs);
    // Return 'true' if the specified 'lhs' and 'rhs' objects do not have the
    // same value and 'false' otherwise.  Two 'Line' objects do not have the
    // same value if any the corresponding values of their 'sequence' and 'name'
    // attributes do not have the same value.

std::istream& operator>>(std::istream& stream, Line& rhs);
    // Assign to the specified 'rhs' object the value extracted from the
    // specified 'stream', and return a reference providing modifiable access
    // to 'stream'.

std::ostream& operator<<(std::ostream& stream, const Line& rhs);
    // Output the value of the specified 'rhs' object to the specified
    // 'stream', and return a reference providing modifiable access to
    // 'stream'.

// ============================================================================
//                        INLINE FUNCTION DEFINITIONS
// ============================================================================

                                // ------------
                                // class Line
                                // ------------

// CREATORS
inline Line::Line()
: orfid()
, evalue()
, line()
{
}

inline
Line::Line(const std::string& name, const std::string& sequence)
: orfid(name)
, evalue(evalue_extractor_from_blast(sequence))
, line(sequence)
{
}

inline
Line::Line(const Line& original)
: orfid(original.orfid)
, evalue(original.evalue)
, line(original.line)
{
}

// MANIPULATORS
inline
Line& Line::operator=(const Line& rhs)
{
    orfid  = rhs.orfid;
    evalue = rhs.evalue;
    line   = rhs.line;

    return *this;
}


inline
void Line::setOrfId(const std::string& value)
{
    orfid.assign(value.begin(), value.end());
}

inline
void Line::setLine(const std::string& value)
{
    line.assign(value.begin(), value.end());
}

inline
void Line::setEvalue(const double& value)
{
	evalue = value;
}

// ACCESSORS
inline const std::string& Line::getOrfId() const
{
    return orfid;
}

inline const std::string& Line::getLine() const
{
    return line;
}


// FREE OPERATORS
inline bool operator==(const Line& lhs, const Line& rhs)
{
    return lhs.getOrfId()  == rhs.getOrfId() && lhs.getLine() == rhs.getLine();
}

inline bool operator!=(const Line& lhs, const Line& rhs)
{
    return lhs.getOrfId() != rhs.getOrfId() || lhs.getLine() != rhs.getLine();
}

// close package namespace

#endif
