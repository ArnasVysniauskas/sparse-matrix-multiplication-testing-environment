#include <assert.h>
#include "application/common/common.h"

struct CompareDoubleTestCase{
    double a;
    double b;
    bool match;
};

void test_compare_double(){
    int no_test_cases = 3;
    struct CompareDoubleTestCase test_cases[] = {
        {0.0, 0.0, true},
        {0.0, 1.0, false},
        {0.0, 111.1 - 111.10, true}
    };

    for (int i = 0; i < no_test_cases; i++)
    {
        assert(
            compare_double(test_cases[i].a, test_cases[i].b) == test_cases[i].match
        );
    }
}

int main(){
    test_compare_double();
}