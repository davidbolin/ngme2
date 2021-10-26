#ifndef NGME_BLOCK_H
#define NGME_BLOCK_H

#include<string>
#include<vector>

#include"latent.h"

using std::string;

class Block
{
friend class Latent;

private:
    string noise;
    
    // std::vector<Latent> latents;
    std::vector<MatrixXd> As, Ks;
public:
    Block(std::vector<Latent>) {}
};




#endif