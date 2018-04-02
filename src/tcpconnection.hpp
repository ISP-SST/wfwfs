#ifndef WFWFS_NETWORK_TCPCONNECTION_HPP
#define WFWFS_NETWORK_TCPCONNECTION_HPP
/***************************************************************************
 *                                                                         *
 *   WFWFS (Wide-field wavefront sensor) software                          *
 *                                                                         *
 *   Copyright (C) 2018-2019 by Tomas Hillberg  <hillberg@astro.su.se>     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 ***************************************************************************/

#include <iostream>
#include <memory>
#include <mutex>
#include <typeinfo>
#include <vector>

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/bind/protect.hpp>

namespace ba = boost::asio;
namespace bip = boost::asio::ip;

namespace wfwfs {

    class TcpConnection : public std::enable_shared_from_this<TcpConnection> {

        static void writeCallback( size_t sent, const boost::system::error_code& error, size_t transferred ) {
            if( !error ) {
                if( sent != transferred ) {
                    std::cerr << "TcpConnection::write: only " << transferred << "/" << sent << " bytes were successfully transferred." << std::endl;
                }
            }
        }

        
    public:

        template <typename T>
        void asyncWrite( const std::shared_ptr<T>& data, size_t sz ) {
            if( mySocket.is_open() ) {
                ba::async_write( mySocket, ba::buffer( data.get(), sz ),
                                            boost::bind( &TcpConnection::writeCallback, sz,
                                                        ba::placeholders::error,
                                                        ba::placeholders::bytes_transferred )
                                        );

            }
        }

        template <typename T>
        void syncWrite( const T* data, size_t sz ) {
            if( mySocket.is_open() ) {
                ba::write( mySocket, ba::buffer( data, sz ) );
            }
        }

        typedef std::shared_ptr<TcpConnection> Ptr;
        typedef std::function<void(Ptr)> callback;

        ~TcpConnection( void );


        static Ptr newPtr( ba::io_service& io_service ) {
            return Ptr( new TcpConnection( io_service ) );
        }

        size_t readline( std::string& line );
        void writeline( const std::string& line );
        template <class T>
        void asyncWrite( const std::vector<T>& data ) {
            size_t sz = data.size();
            if (!sz) return;
            std::shared_ptr<T> tmp( new T[sz], [](T* p){ delete[] p;} );
            memcpy(tmp.get(),data.data(),sz*sizeof(T));
            asyncWrite(tmp, sz*sizeof(T));
        }

        template <class T>
        void asyncWrite( const T& data ) {
            asyncWrite(std::make_shared<T>(data), sizeof(T) );
        }

        template <class T>
        void syncWrite( const std::vector<T>& in ) {
            size_t sz = in.size();
            if (!sz) return;
            syncWrite( in.data(), sz*sizeof(T) );
        }

        template <class T>
        void syncWrite( const T& data ) {
            syncWrite( &data, sizeof(T) );
        }

        bip::tcp::socket& socket() { return mySocket; }
        operator bool() const { return mySocket.is_open(); };

        void connect( std::string host, std::string service );
        void close( void );
        callback getCallback( void ) { std::unique_lock<std::mutex> lock(mtx); return activityCallback; };
        bool hasUrgentCallback( void ) const { return (urgentCallback != nullptr); };
        void setCallback( callback cb = nullptr ) { std::unique_lock<std::mutex> lock(mtx); activityCallback = cb; };
        void setUrgentCallback( callback cb = nullptr ) { std::unique_lock<std::mutex> lock(mtx); urgentCallback = cb; };
        void setErrorCallback( callback cb = nullptr ) { std::unique_lock<std::mutex> lock(mtx); errorCallback = cb; };
        void uIdle( void );
        void urgentHandler( const boost::system::error_code& error, size_t transferred );
        void idle( void );
        void onActivity( const boost::system::error_code& error );
        void setSwapEndian(bool se) { swapEndian_ = se; };
        bool getSwapEndian(void) { return swapEndian_; };
        
        unsigned long getRemoteIP( void );
        unsigned short getRemotePort( void );
        
        void lock(void) { mtx.lock(); };
        void unlock(void) { mtx.unlock(); };
        bool try_lock(void) { return mtx.try_lock(); };
        
        void sendUrgent( uint8_t c );
        void receiveUrgent( uint8_t& c );
        uint8_t getUrgentData(void) { return urgentData; };
        
        TcpConnection& operator<<( const uint8_t& );
        TcpConnection& operator>>( uint8_t& );
        
        TcpConnection& operator<<( const std::vector<std::string>& );
        TcpConnection& operator>>( std::vector<std::string>& );

        template <typename T>
        TcpConnection& operator<<( const T& in ) {
            syncWrite(&in, sizeof(T));
            return *this;
        }


        template <typename T>
        TcpConnection& operator>>( T& out ) {
            if( ba::read( mySocket, ba::buffer( &out, sizeof(T) ) ) < sizeof(T) ) {
                out = T();
                throw std::ios_base::failure( std::string("TcpConnection: Failed to receive ")+typeid(T).name() );
            }
            return *this;
        }

    private:
        TcpConnection( ba::io_service& io_service );
        TcpConnection( const TcpConnection& ) = delete;

        callback activityCallback;
        callback urgentCallback;
        callback errorCallback;
        bip::tcp::socket mySocket;
        ba::io_service& myService;
        bool swapEndian_;
        uint8_t urgentData;
        bool urgentActive;
        std::mutex mtx;
        std::vector<std::string> lines;

    public:

    };

}   // wfwfs

#endif // WFWFS_NETWORK_TCPCONNECTION_HPP
