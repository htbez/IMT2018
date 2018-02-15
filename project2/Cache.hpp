#ifndef cache_hpp
#define cache_hpp
#include <functional>
#include <map>

template<class T1, class T2>
class Cache
{
    public :    
        Cache () {}
        Cache ( const std::function<T2(T1)> f_ ){
            f = f_ ;
        }        
        T2 operator() (T1 key ) const {
            auto index = memory_map.find(key) ;
            if ( index != memory_map.end() )
                return memory_map[key] ;
            else
                return calculateData( key ) ;
        }        
        
        void setf(const std::function<T2(T1)> f_ ) {
            f = f_;
        }
        
        void clear() {
            memory_map.clear ;
        }        
        void erase(T1 key ){
            memory_map.erase( key ) ;
        }
        
    private:
        std::function<T2(T1)> f ;        
        mutable std::map<T1,T2> memory_map ;
        
        T2 calculateData (T1 key ) const {
            put ( key , f( key ) ) ;
            return memory_map[key] ;
        }
        
        template<class T3>
        void put (T1 key , T3 value ) const {
            memory_map[key] = static_cast<T2>(value) ;
        }
};
#endif
