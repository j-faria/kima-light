#ifndef DNest4_Prior
#define DNest4_Prior

#include <type_traits>

template <class T>
class prior {
  public:
    T dist;

    // template<typename... args> 
    // prior(args... a) : dist(a...) {}; 

    prior() : dist() {};
    prior(double arg) : dist(arg) {};
    prior(double arg1, double arg2) : dist(arg1, arg2) {};

    // class D { public: operator C() { return c; }  C c; };
    // template <class F>
    // operator F() { return this; };

    // template <class F>
    // F& operator= (F f){
      // F newF = f;
      // std::remove_cv<f>::type type1
      // return f;
    // }

    // template <typename U>
    // prior(std::initializer_list il): dist(il) {}
    // prior(std::initializer_list<U> il): v(il) {}
    // prior(double ... args): dist(&args...) {}
};

#endif

