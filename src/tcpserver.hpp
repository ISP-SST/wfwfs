#ifndef WFWFS_NETWORK_TCPSERVER_HPP
#define WFWFS_NETWORK_TCPSERVER_HPP
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

namespace wfwfs {

    class TcpServer : private boost::noncopyable {

    public:

        TcpServer( boost::asio::io_service& io_service, uint16_t port );

        void accept(void);
        void setCallback( TcpConnection::callback cb = nullptr ) { onConnected = cb; };

    private:

        void onAccept( TcpConnection::Ptr conn, const boost::system::error_code& error );

        bip::tcp::acceptor acceptor;
        TcpConnection::callback onConnected;

    };

}   // wfwfs

#endif // WFWFS_NETWORK_TCPSERVER_HPP
