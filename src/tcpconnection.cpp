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
#include "tcpconnection.hpp"

#include "datautil.hpp"
#include "stringutil.hpp"

#include <thread>

#include <boost/algorithm/string.hpp>

namespace ba = boost::asio;

using namespace wfwfs;
using namespace std;

#ifdef DEBUG_
//#define DBG_NET_
#endif


namespace {

#ifdef DBG_NET_
    static atomic<int> connCounter(0);
#endif

}


TcpConnection::TcpConnection( boost::asio::io_service& io_service )
    : activityCallback( nullptr ), urgentCallback( nullptr ), errorCallback( nullptr ), mySocket( io_service ),
    myService( io_service ), swapEndian_(false), urgentActive(false) {
#ifdef DBG_NET_
    cout << "Constructing TcpConnection: (" << hexString(this) << ") new instance count = " << (connCounter.fetch_add(1)+1) << endl;
#endif
}

            
TcpConnection::~TcpConnection( void ) {
    mySocket.close();
#ifdef DBG_NET_
    cout << "Destructing TcpConnection: (" << hexString(this) << ") new instance count = " << (connCounter.fetch_sub(1)-1) << endl;
#endif
}


namespace {
    const string delimiter = "\r\n";
}

size_t TcpConnection::readline( string& out ) {

    out.clear();
    
    lock_guard<mutex> lock(mtx);

    boost::asio::socket_base::bytes_readable command(true);
    mySocket.io_control(command);
    size_t bytes_readable = command.get();
    if( bytes_readable ) {
        try {
            boost::asio::streambuf streambuf;
            boost::asio::read_until( mySocket, streambuf, delimiter );
            string line{ buffers_begin(streambuf.data()),
                         buffers_begin(streambuf.data()) + streambuf.size() - delimiter.size()};
            streambuf.consume( streambuf.size() );
            if( line.size() ) {
                vector<string> split_line;
                boost::iter_split( split_line, line, boost::first_finder(delimiter) );
                for( const auto& l: split_line ) {
                    if( !l.empty() ) {
                        lines.push_back( l );
                    }
                }
                
            }
            
        } catch( const exception& e ) {
#ifdef DBG_NET_
            cerr << "TcpConnection::readline  e = " << e.what() << endl
#endif
        }
    }

    while( lines.size() && out.empty() ) {
        out = lines.front();
        lines.erase( lines.begin() );
    }

    return out.size();

}


void TcpConnection::writeline( const string& line ) {

    boost::asio::write( socket(), boost::asio::buffer(line + delimiter) );

}


void TcpConnection::connect( string host, string service ) {
    
    if( host == "" ) host = "localhost";

    try {
        ba::ip::tcp::resolver::query query( host, service );
        ba::ip::tcp::resolver resolver( myService );
        ba::ip::tcp::resolver::iterator destination = resolver.resolve( query );
        ba::ip::tcp::resolver::iterator end ;
        ba::ip::tcp::endpoint endpoint;

        while( destination != end ) {
            try {
                mySocket.connect( *destination++ );
            } catch ( const boost::system::system_error& ) {
                mySocket.close();
            }
            if( mySocket.is_open() ) return;
        }
    } catch ( ... ) {
        // TODO 
    }
    mySocket.close();

}


void TcpConnection::close( void ) {

    unique_lock<mutex> lock(mtx);
    activityCallback = nullptr;
    urgentCallback = nullptr;
    errorCallback = nullptr;
    mySocket.close();
    
}


void TcpConnection::uIdle( void ) {

    unique_lock<mutex> lock(mtx);
    if( !urgentCallback || urgentActive ) {
        return;
    }

    urgentActive = true;
    
    if( mySocket.is_open() ) {
        lock.unlock();
                
        mySocket.async_receive( boost::asio::buffer( &urgentData, 1 ), ba::socket_base::message_out_of_band,
                                  boost::bind( &TcpConnection::urgentHandler, this,
                                               ba::placeholders::error, ba::placeholders::bytes_transferred) );
    }

}


void TcpConnection::urgentHandler( const boost::system::error_code& error, size_t transferred ) {
    
    try {
        unique_lock<mutex> lock(mtx);
        if( !urgentActive ) return;     // prevent multiple uIdle
        urgentActive = false;
        if( mySocket.is_open() && !mySocket.at_mark() ) {
            lock.unlock();
            uIdle();
            return;
        }
    } catch(...) {
        return;
    }
    
    if( !error ) {
        if( urgentCallback ) {
            if( mySocket.is_open() ) {
                std::thread( urgentCallback, shared_from_this() ).detach();
            }
        }
    } else {
        if( ( error == ba::error::eof ) || ( error == ba::error::connection_reset ) ) {
            mySocket.close();
        } else {
            if( errorCallback ) {
                std::thread( errorCallback, shared_from_this() ).detach();
            } else {
                //throw std::ios_base::failure( "TcpConnection::urgentHandler: error: " + error.message() );
            }
        }
    }

}

void TcpConnection::idle( void ) {

    unique_lock<mutex> lock(mtx);
    if( !activityCallback ) {
        return;
    }

    if( mySocket.is_open() ) {
        lock.unlock();
        mySocket.async_read_some( ba::null_buffers(),
                                  boost::bind( &TcpConnection::onActivity, this, ba::placeholders::error ) );

    }
    
}

void TcpConnection::onActivity( const boost::system::error_code& error ) {

    unique_lock<mutex> lock(mtx);
    if( !error ) {
        if( activityCallback ) {
            if( mySocket.is_open() ) {
                std::thread( activityCallback, shared_from_this() ).detach();
            }
        }
    } else {
        if( ( error == ba::error::eof ) || ( error == ba::error::connection_reset ) ) {
            mySocket.close();
        } else {
            if( errorCallback ) {
                std::thread( errorCallback, shared_from_this() ).detach();
            } else {
                throw std::ios_base::failure( "TcpConnection::onActivity: error: " + error.message() );
            }
        }
    }

}

unsigned long TcpConnection::getRemoteIP( void ) {
    unsigned long ip(0);
    try {
        bip::address address = mySocket.remote_endpoint().address();
        if( address.is_v4() ) {
            ip = address.to_v4().to_ulong();
        } else if( address.is_v6() ) {
            ip = address.to_v4().to_ulong();    // TODO fix v6 version
        }
    } catch ( const boost::system::system_error& e ) {
        mySocket.close();
    }
    return ip;
}


unsigned short TcpConnection::getRemotePort( void ) {
    unsigned short port(0);
    try {
        port = mySocket.remote_endpoint().port();
    } catch ( const boost::system::system_error& e ) {
        mySocket.close();
    }
    return port;
}


void TcpConnection::sendUrgent( uint8_t c ) {
    if( mySocket.is_open() ) {
        mySocket.send( boost::asio::buffer( &c, 1 ), ba::socket_base::message_out_of_band );
    }
}

 
void TcpConnection::receiveUrgent( uint8_t& c ) {
    if( mySocket.is_open() && mySocket.at_mark() ) {
        mySocket.receive( boost::asio::buffer( &c, 1 ), ba::socket_base::message_out_of_band );
    }
}

 
TcpConnection& TcpConnection::operator<<( const uint8_t& in ) {
    syncWrite(&in, sizeof(uint8_t));
    return *this;
}


TcpConnection& TcpConnection::operator>>( uint8_t& out ) {
    if( ba::read( mySocket, ba::buffer( &out, sizeof( uint8_t ) ) ) < sizeof( uint8_t ) ) {
        out = 0;
        throw std::ios_base::failure( "Failed to receive command." );
    }
    return *this;
}


TcpConnection& TcpConnection::operator<<( const std::vector<std::string>& in ) {
    uint64_t inSize(0);
    if( in.size() ) {
        uint64_t messagesSize = wfwfs::size( in );
        size_t totSize = messagesSize+sizeof(uint64_t);
        shared_ptr<char> tmp( new char[totSize], []( char* p ){ delete[] p; } );
        char* ptr = tmp.get();
        ptr += wfwfs::pack( ptr, messagesSize );
        wfwfs::pack( ptr, in );
        syncWrite(tmp.get(), totSize);
    } else {
        syncWrite( inSize );
    }
    return *this;
}


TcpConnection& TcpConnection::operator>>( std::vector<std::string>& out ) {
    uint64_t blockSize, received;
    received = boost::asio::read( mySocket, boost::asio::buffer( &blockSize, sizeof(uint64_t) ) );
    if( received == sizeof(uint64_t) ) {
        if( swapEndian_ ) swapEndian( blockSize );
        if( blockSize ) {
            shared_ptr<char> buf( new char[blockSize+1], []( char* p ){ delete[] p; } );
            char* ptr = buf.get();
            memset( ptr, 0, blockSize+1 );
            received = boost::asio::read( mySocket, boost::asio::buffer( ptr, blockSize ) );
            if( received ) {
                unpack( ptr, out, swapEndian_ );
            }
        }
    }
    return *this;
}

