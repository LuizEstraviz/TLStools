#ifndef TLSTOOLSBASE_H
#define TLSTOOLSBASE_H
#include <string>

using namespace std;

class TLStoolsBase
{
    public:
        /** Default constructor */
        TLStoolsBase();
        /** Default destructor */
        virtual ~TLStoolsBase();
        virtual string getHelp() = 0;

    protected:

    private:
};

#endif // TLSTOOLSBASE_H
