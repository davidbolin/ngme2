#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

enum Color { red, green, blue };

int main(int argc, char const *argv[])
{
    vector<bool> v (10, true);
    v.at(8) = false;
    bool all_true = false;

    if (std::find(begin(v), end(v), false) == end(v))
        all_true = true;

    cout << green << endl;
}