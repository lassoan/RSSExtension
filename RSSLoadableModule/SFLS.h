#ifndef SFLS_h_
#define SFLS_h_


// std
#include <list>

template<typename T>
class SFLSNode
{
    public:
    SFLSNode(T x, T y, T z) {SFLSNodeComponent1 = x; SFLSNodeComponent2 = y; SFLSNodeComponent3 = z;}

    T SFLSNodeComponent1;
    T SFLSNodeComponent2;
    T SFLSNodeComponent3;
};

class CSFLS
{
public:
    typedef CSFLS Self;

    typedef SFLSNode<long> NodeType;
    typedef std::list< NodeType > CSFLSLayer;

    // ctor
    CSFLS() {}

    CSFLSLayer m_lz;
    CSFLSLayer m_ln1;
    CSFLSLayer m_ln2;
    CSFLSLayer m_lp1;
    CSFLSLayer m_lp2;
};

#endif
