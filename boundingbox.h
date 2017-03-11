#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include <array> // std::array
#include <vector> // std::array
//#include <CL/cl_platform.h>
#include <limits> // numeric_limits
#include <ostream>

template<typename Container, typename T>
class ConstIteratorProxy
{
public:
    ConstIteratorProxy(const Container& container) : container_(container) {}
    typename Container::const_iterator begin() const {
        return container_.begin();
    }
    typename Container::const_iterator end() const {
        return container_.end();
    }
    T operator[](size_t pos) const {
        return container_.at(pos);
    }
    friend std::ostream& operator<<(std::ostream& stream,
                                    const ConstIteratorProxy& proxy) {
        for (auto& element : proxy.container_)
            stream << element << " ";
        return stream;
    }

private:
    const Container& container_;
};

template<typename T, unsigned int N>
class BoundingBox
{
private:
    std::array<T, N> minVertex_;
    std::array<T, N> maxVertex_;
public:
    BoundingBox();
    std::vector<std::array<T, N> > getBoundingBox() const;
    void setMinVertex(T vertice[N]);
    void setMinVertex(std::array<T, N>& vertice);
    void setMaxVertex(T vertice[N]);
    void setMaxVertex(std::array<T, N>& vertice);

    //Sets minimum bounding box
    void setBoundingBox(std::vector<std::array<T, N> > &vertices);
#ifdef __CL_PLATFORM_H
    void setBoundingBox(std::vector<cl_float4> &vertices);
#endif

    //void setInnerBoundingBox(std::vector<cl_float4>);

    ConstIteratorProxy<std::array<T, N>, T> minVertex() const {
        return ConstIteratorProxy<std::array<T, N>, T>(minVertex_);
    }
    ConstIteratorProxy<std::array<T, N>, T> maxVertex() const {
        return ConstIteratorProxy<std::array<T, N>, T>(maxVertex_);
    }
};


template<typename T, unsigned int N>
BoundingBox<T, N>::BoundingBox()
{
    minVertex_.fill(-std::numeric_limits<T>::max());
    maxVertex_.fill(std::numeric_limits<T>::max());
}

template<typename T, unsigned int N>
std::vector<std::array<T, N> > BoundingBox<T, N>::getBoundingBox() const
{
    std::vector<std::array<T, N> > box {minVertex_, maxVertex_};
    return box;
}

template<typename T, unsigned int N>
void BoundingBox<T, N>::setMinVertex(T vertex[N])
{
    for (unsigned int i = 0; i < N; i++)
        minVertex_.at(i) = vertex[i];
}

template<typename T, unsigned int N>
void BoundingBox<T, N>::setMinVertex(std::array<T, N>& vertex)
{
    for (unsigned int i = 0; i < N; i++)
        minVertex_.at(i) = vertex[i];
}

template<typename T, unsigned int N>
void BoundingBox<T, N>::setMaxVertex(T vertex[N])
{
    for (unsigned int i = 0; i < N; i++)
        maxVertex_.at(i) = vertex[i];
}

template<typename T, unsigned int N>
void BoundingBox<T, N>::setMaxVertex(std::array<T, N>& vertex)
{
    for (unsigned int i = 0; i < N; i++)
        maxVertex_.at(i) = vertex[i];
}

template<typename T, unsigned int N>
void BoundingBox<T, N>::setBoundingBox(std::vector<std::array<T, N>>&
                                            vertices)
{
    T min[N];
    for (unsigned int i = 0; i < N; i++)
        min[i] = std::numeric_limits<T>::max();

    T max[N];
    for (unsigned int i = 0; i < N; i++)
        max[i] = -std::numeric_limits<T>::max();

    for (auto& vertex : vertices) {
        for (unsigned int i = 0; i < N; i++) {
            if (min[i] > vertex[i]) min[i] = vertex[i];
            if (max[i] < vertex[i]) max[i] = vertex[i];
        }
    }
    setMinVertex(min);
    setMaxVertex(max);
}

#ifdef __CL_PLATFORM_H
template<typename T, unsigned int N>
void BoundingBox<T, N>::setBoundingBox(std::vector<cl_float4> &vertices)
{
    T min[N];
    for (unsigned int i = 0; i < N; i++)
        min[i] = std::numeric_limits<T>::max();

    T max[N];
    for (unsigned int i = 0; i < N; i++)
        max[i] = -std::numeric_limits<T>::max();

    for (auto& vertex : vertices) {
        for (unsigned int i = 0; i < N; i++) {
            if (min[i] > vertex.s[i]) min[i] = vertex.s[i];
            if (max[i] < vertex.s[i]) max[i] = vertex.s[i];
        }
    }
    setMinVertex(min);
    setMaxVertex(max);
}

#endif


#endif // BOUNDINGBOX_H
