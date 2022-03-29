#pragma once

#include "db/Database.h"
#include "gr_db/GrDatabase.h"
#include "single_net/SingleNetRouter.h"
#include "db/RsynService.h"
#include "db/Cell.h"
#include "db/GeoPrimitive.h"

namespace crp{

struct Node{
        int idx,r,s;
        double cost;
    };
struct Edge{
    int i , j;
};

class Solver{
public:
    Solver();
    ~Solver();

    void weightBuilder();
    

};//end Solver
};//end namespace crp
