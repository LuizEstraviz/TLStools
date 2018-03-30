#include <iostream>

using namespace std;

int main()
{

 #ifdef DEBUGMODE
    cout << "dbg..." << endl;
 #endif // DEBUGMODE

    cout << "Hello world!" << endl;
    return 0;
}
