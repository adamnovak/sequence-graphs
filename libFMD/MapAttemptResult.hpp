#ifndef MAPATTEMPTRESULT_HPP
#define MAPATTEMPTRESULT_HPP

/**
 * A triple to hold the return values from FMD::mapPosition() or
 * FMD::partialMap(). Holds a flag for whether the mapping succeeded or not, an
 * FMDPosition corresponding either to where the character mapped or the longest
 * search starting at the character that did actually return results, and the
 * number of characters in the FMDPosition's search pattern.
 */
struct MapAttemptResult
{
    bool is_mapped;
    FMDPosition position;
    size_t characters;
};

struct creditMapAttemptResult
{
    bool is_mapped;
    FMDPosition position;
    size_t characters;
    size_t maxCharacters;
};

#endif
