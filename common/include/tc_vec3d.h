#pragma once

#include <math.h>
#include "adolc/adolc.h"

template<typename a, typename b> struct select_type{ typedef adouble type; };
template<> struct select_type<double, double>{ typedef double type; };

#define tmpl_T template<typename T>
#define tmpl_T1_T2 template<typename T1, typename T2>
#define tmpl_type typename select_type<T1, T2>::type
#define FOR(ii,n) for(int ii=0; ii<n; ii++)

tmpl_T void set_vec3d(T* v1, const T* v2)
{ FOR(i,3) v1[i]=v2[i]; }

tmpl_T void add_vec3d(T* v1, const T* v2)
{ FOR(i,3) v1[i]+=v2[i]; }

tmpl_T void subtract_vec3d(T* v1, const T* v2)
{ FOR(i,3) v1[i]-=v2[i]; }

tmpl_T void divide_vec3d(T* v1, const T a)
{ FOR(i,3) v1[i]/=a; }

tmpl_T T vec3dNorm(T* v)
{ return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

tmpl_T T normVec3d(T* v)
{ T m=vec3dNorm(v); FOR(i,3) v[i]=v[i]/m; return m; }

tmpl_T T normVec3d(T* vN, const T* v)
{ T m=vec3dNorm(v); FOR(i,3) vN[i]=v[i]/m; return m; }

tmpl_T void normalize_vec3d(T* v)
{ T m=normVec3d(v); }

tmpl_T void vecMinVec3d(T *w, const T *u, const T *v)
{ FOR(i,3) w[i]=u[i]-v[i]; }

tmpl_T1_T2 tmpl_type vecDotVec3d(const T1* v1, const T2* v2)
{ tmpl_type d; FOR(i,3) d+=v1[i]*v2[i]; return d; }

tmpl_T void matMultMatR5(T (*w)[5], T (*u)[5], T (*v)[5])
{ FOR(i,5) FOR(j,5) { w[i][j]=0; FOR(k,5) w[i][j]+=u[i][k]*v[k][j]; } }

tmpl_T void matMultDiagMatR5(T (*w)[5], T (*u)[5], T *v)
{ FOR(i,5) FOR(j,5) w[i][j]=u[i][j]*v[j]; }

tmpl_T void matMultVecR5(T *w, T (*u)[5], T *v)
{ FOR(i,5) { w[i]=0; FOR(j,5) w[i]+=u[i][j]*v[j]; } }

#undef tmpl_T
#undef tmpl_T1_T2
#undef tmpl_type
#undef FOR
