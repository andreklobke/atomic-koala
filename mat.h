//    Created by andre klobke on 24.11.13.
//    Copyright (c) 2013 andre klobke. All rights reserved.
//
//    "atomic-koala" is a free library: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This code is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this file.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __ak__mat__
#define __ak__mat__

#include <sstream>
#include <string>
#include <ostream>
#include <array>
#include <numeric>          // std::inner_product
#include <functional>       // std::plus
#include <algorithm>        // std::transform, std::fill
#include <math.h>           // sqrt, pow
#include <type_traits>      // std::is_integral
#include <initializer_list> // std::initializer_list

namespace ak{
    template<size_t ROW=4, size_t COL=ROW,typename T=float> class mat;
    template<size_t N=3, typename T=float> class vec;
    
//typedefs
    typedef mat<2, 2, float> mat2;
    typedef mat<3, 3, float> mat3;
    typedef mat<4, 4, float> mat4;
    typedef mat<5, 5, float> mat5;
    typedef vec<2, float> vec2;
    typedef vec<3, float> vec3;
    typedef vec<4, float> vec4;
    typedef vec<5, float> vec5;
    
//common_matrix_type_operations
    namespace _common_matrix_type_operations{
        struct TAG{};
        template<typename M>M operator+(const M& left, const M& right){return M(left)+=right;}
        template<typename M>M operator-(const M& left, const M& right){return M(left)-=right;}
        template<typename M>M operator*(const M& left, const M& right){return M(left)*=right;}
        template<typename M, typename T>M operator+(const M& left, const T& right){return M(left)+=right;}
        template<typename M, typename T>M operator-(const M& left, const T& right){return M(left)-=right;}
        template<typename M, typename T>M operator*(const M& left, const T& right){return M(left)*=right;}
    }

//mat
    template<size_t ROW, size_t COL,typename T>
    class mat :_common_matrix_type_operations::TAG{
    private:
        template<size_t H, size_t W, typename D>
        struct static_constructor{};
        template<size_t W, typename D>
        struct static_constructor<W,W,D>{
            mat<W,W,D> ID;
            static_constructor(void){
                D ar[W*W] = {static_cast<T>(0)};
                for(int i=0; i<W*W; i+=W+1)
                    ar[i]= static_cast<T>(1);
                ID = mat<W,W,D>(ar); 
            }
        };
    protected:
        std::array<T, ROW*COL> data;
        
        template<typename D>
        void init(const D (&data) [ROW*COL]){
            for(int i=0; i<ROW*COL; ++i)
                this->data[i] = static_cast<T>(data[i]);
        }
    public:
//        static const size_t HEIGHT = ROW;
//        static const size_t WIDTH = COL;
        static static_constructor<ROW,COL,T> CONST;
        bool my;
        inline T* begin(){return data.begin();}
        inline T* end(){return data.end();}
        inline const T* cbegin()const{return data.cbegin();}
        inline const T* cend()const{return data.cend();}
        inline size_t size() const {return ROW*COL;}
       
        mat(){}
        mat(const T&& arg){
            std::fill(begin(),end(), static_cast<T>(arg));
        }
        mat(const std::initializer_list<T> args){
            assert(args.size()==ROW*COL && "wrong number of arguments");
            std::copy(args.begin(), args.end(), data.begin());
        }
        template<typename D>
        mat(const D (&data) [ROW*COL]){
            init<D>(data);
        }
        inline T at(const size_t row, const size_t col=0) const{
            return data[row*COL+col];
        }
        inline T& at(const size_t row, const size_t col=0){
            return data[row*COL+col];
        }
        std::string toString() const{
            std::stringstream ss;
            for(int i=0; i<ROW; i++){
                ss<<(i==0?'[':' ');
                for(int j=0; j<COL; j++){
                    ss<<this->data[i*COL+j];
                    ss<<(i==ROW-1 && j==COL-1 ? ']' : ',');
                }
                ss<<'\n';
            }
            return ss.str();
        }
        mat<COL,ROW,T> transpose(){
            T res[ROW*COL];
            for(int i=0; i<ROW; i++)
                for(int j=0; j<COL; j++)
                    res[j*ROW+i] =this->data[i*COL+j];
            return ak::mat<COL,ROW,T>(res);
        }
        bool operator==(const ak::mat<ROW,COL,T>& right){
            if(this==&right)
                return true;
            for(int i=0; i<ROW*COL; i++)
                if(data[i]!=right.data[i])
                    return false;
            return true;
        }
        bool operator!=(const ak::mat<ROW,COL,T>& right){
            return !(*this == right);
        }
        template<size_t COL_RIGHT>
        mat<ROW,COL_RIGHT,T> operator*(const mat<COL,COL_RIGHT,T>& right){
            T res[ROW*COL_RIGHT];
            for(int i=0; i<ROW; i++)
                for(int k=0; k<COL_RIGHT; k++){
                    res[i*COL_RIGHT+k]=0;
                    for(int j=0; j<COL; j++)
                        res[i*COL_RIGHT+k]+= at(i,j)*right.at(j,k);
                }
            return mat<ROW,COL_RIGHT,T>(res);
        }
        mat<ROW,COL,T> operator+=(const mat<ROW,COL,T>& right){
            std::transform(cbegin(), cend(), right.cbegin(), begin(), std::plus<T>());
            return *this;
        }
        mat<ROW,COL,T> operator-=(const mat<ROW,COL,T>& right){
            std::transform(cbegin(), cend(), right.cbegin(), begin(), std::minus<T>());
            return *this;
        }
        mat<ROW,COL,T> operator*=(const T& right){
            std::for_each(begin(),end(),[&](T& elem){elem *= right;});
            return *this;
        }
        mat<ROW,COL,T> operator+=(const T& right){
            std::for_each(begin(),end(),[&](T& elem){elem += right;});
            return *this;
        }
        mat<ROW,COL,T> operator-=(const T& right){
            std::for_each(begin(),end(),[&](T& elem){elem -= right;});
            return *this;
        }
        mat<ROW,COL,T> operator-(){
            return mat<ROW,COL>(*this)*=static_cast<T>(-1.f);
        }
        //T* operator[](int row){ return data.begin() + row*COL; }
    };
    template<size_t H, size_t W, typename D>
    mat<H,W,D>::static_constructor<H, W, D> mat<H,W,D>::CONST;
    
//vec
    template<size_t N, typename T>
    class vec: public mat<N,1,T>{
    private:
        template<size_t M,typename D>
        struct static_constructor{
            vec<M,D> ZERO,ONE;
            static_constructor(void):
            ZERO(vec<M,D>(static_cast<D>(0))),
            ONE(vec<M,D>(static_cast<D>(1))) {}
        };
        template<typename D>
        struct static_constructor<3,D> {
            vec<3,D> X,Y,Z,NEG_Z,ZERO,ONE;
            static_constructor(void):
            X(vec<3,D>{1,0,0}),
            Y(vec<3,D>{0,1,0}),
            Z(vec<3,D>{0,0,1}),
            NEG_Z(vec<3,D>{0,0,-1}),
            ZERO(vec<3,D>(static_cast<D>(0))),
            ONE(vec<3,D>(static_cast<D>(1))){}
        };
    public:
        static static_constructor<N,T> CONST;
        
        vec(){}
        vec(const T&& arg) : mat<N,1,T>(std::move(arg)){}
        //vec(int arg) : mat<N,1,T>(std::move(arg)){}
        template<typename D> vec(const D (&data) [N]) : mat<N,1,T>(data){}
        vec(const mat<N,1,T>& a0) :mat<N,1,T>(a0){}
        vec(const vec<N-1,T>& a0, T a1) : mat<N,1,T>(){
            std::copy(a0.begin(), a0.end(), this->begin());
            (this->at(N-1))=a1;
        }
        vec(T a0, const vec<N-1,T>& a1) : mat<N,1,T>(){
            std::copy(a1.begin(), a1.end(), this->begin()+1);
            (this->at(0))=a0;
        }
        vec(const std::initializer_list<T> args) : mat<N,1,T>(args){ }

        inline const T len() const {
            //return sqrt( std::accumulate(this->cbegin(),this->cend(),static_cast<T>(0),[](T acc, T n){return acc+n*n;}) );
            return sqrt(std::inner_product(this->cbegin(), this->cend(), this->cbegin(), static_cast<T>(0)));
        }
        vec<N,T> normalize() const{
            T len = static_cast<T>(1)/this->len();
            return vec<N,T>(*this) *= len;
        }
        vec<N,T> operator*=(const mat<N,1,T>& right){
            std::transform (this->cbegin(), this->cend(), right.cbegin(), this->begin(), std::multiplies<T>());
            return *this;
        }
        vec<N,T> operator+=(const mat<N,1,T>& right){mat<N,1,T>::operator+=(right);return *this;}
        vec<N,T> operator-=(const mat<N,1,T>& right){mat<N,1,T>::operator-=(right);return *this;}
        vec<N,T> operator*=(const T& right){mat<N,1,T>::operator*=(right);return *this;}
        vec<N,T> operator+=(const T& right){mat<N,1,T>::operator+=(right);return *this;}
        vec<N,T> operator-=(const T& right){mat<N,1,T>::operator-=(right);return *this;}
        vec<N,T> operator-(){return mat<N,1,T>::operator-();}
        //T* operator[](int n){ return this->data[n]; }
   
    };//eof vec
    template<size_t N, typename T>
    vec<N,T>::static_constructor<N,T> vec<N,T>::CONST;
    
//convenience methods
    template<size_t ROW, size_t COL, typename T> typename std::enable_if<ROW==1||COL==1,
    T>::type dot(const mat<ROW,COL,T>& left, const mat<ROW,COL,T>& right){
        return std::inner_product(left.cbegin(), left.cend(), right.cbegin(), static_cast<T>(0));
    }
    template<typename T>
    vec<3,T> cross(const ak::mat<3,1,T>& left,const ak::mat<3,1,T>& right){
        T res[3];
        res[0] = (left.at(1)*right.at(2)) - (left.at(2)*right.at(1));
        res[1] = (left.at(2)*right.at(0)) - (left.at(0)*right.at(2));
        res[2] = (left.at(0)*right.at(1)) - (left.at(1)*right.at(0));
        return vec<3,T>(res);
    }
    template<size_t ROW, size_t COL, typename T>
    std::ostream& operator<<( std::ostream& stream, const mat<ROW,COL,T> mat){
        for(int i=0; i<ROW; i++){
            stream<<(i==0?'[':' ');
            for(int j=0; j<COL; j++){
                if(j==0)stream<<'{';
                stream<<mat.at(i,j);//mat.data[i*COL+j];
                if(j==COL-1)stream<<'}';
                stream<<(i==ROW-1 && j==COL-1 ? ']' : ',');
            }
        }
        return stream;
    }
    
}//eof ak
#endif