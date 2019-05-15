#ifndef WFWFS_DAEMON_HPP
#define WFWFS_DAEMON_HPP
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

#include "application.hpp"
#include "array.hpp"
#include "camera.hpp"
#include "datautil.hpp"
#include "seeing.hpp"
#include "filefits.hpp"
#include "tcpserver.hpp"
#include "seeing.hpp"


#include <mutex>
#include <thread>
#include <vector>

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>

namespace wfwfs {
    
    struct Host {
        unsigned long address;
        unsigned short port;
        bool operator<( const Host& rhs ) {
            if( address != rhs.address ) return address < rhs.address;
            return port < rhs.port;
        }
    };
    inline std::ostream& operator<<(std::ostream& os, const Host& h) {
        os << std::to_string(h.address) << ":" << std::to_string(h.port);
        return os;
    }


    class Daemon : public wfwfs::Application {

    public:

        explicit Daemon ( boost::program_options::variables_map& vm );
        virtual ~Daemon ( void );
        
        void reset(void);
        void restart(void);
        void stop(void);
        
    private:

        void init(void);
        //void init_cells( const std::vector<PointI>& );
        void parsePropertyTree( boost::property_tree::ptree& );
        
        void queue_frame( const uint8_t*, boost::posix_time::ptime );
        void get_frame( TcpConnection::Ptr conn, int x1, int y1, int x2, int y2, int scale=1, bool darkflat=false, size_t fsel=0, bool do_histo=false );
        std::string make_filename( int, long, const std::string& state );
        void save_burst( int, int, long, std::string );
        template <typename T> void df_cell( Frame& f, Cell& c );
        void copy_cell_data( Frame& f );
        
        void threadLoop( void );
        void maintenance( void );
        bool doWork(void);
        void setThreads( int nThreads );
        
        void connect( Host& host, TcpConnection::Ptr& conn );
        
        void connected( TcpConnection::Ptr );
        void onMessage( TcpConnection::Ptr );
        bool processCmd( TcpConnection::Ptr, const std::string& );
        
        void addConnection(const Host&, TcpConnection::Ptr&);
        void removeConnection(TcpConnection::Ptr);

        void subscribe( const TcpConnection::Ptr& conn, std::string tag );
        void unsubscribe( const TcpConnection::Ptr& conn, std::string tag );
        void broadcast( std::string tag, std::string message, TcpConnection::Ptr skip=TcpConnection::Ptr() );
        
        void start_cam( void );
        void stop_cam( void );

        void play( void );
        void pause( void );
        
        void light( bool );

        std::string list_calib( void );
        std::string list_calibs( void );
        void load_dark( int id=-1 );
        void load_flat( int id=-1 );
        
        void updateCalibID( void );
        void makeHdr( std::shared_ptr<Fits>& );
        int accumulate( Array<uint32_t>&, size_t );
        void do_accumulation( Array<float>&, Array<uint32_t>&, std::string, size_t );
        void darkburst( size_t );
        void flatburst( size_t );
        void update_gain( void );
        
        void die(void);
        void die( TcpConnection::Ptr&, bool urgent=false );
        void softExit(void);
        void reset( TcpConnection::Ptr&, bool urgent=false );
        
        void thumbnail( TcpConnection::Ptr );
        
        int test_threads( int nT, int t=10 );
        
        boost::program_options::variables_map& settings;
        
        uint16_t nQueuedJobs;
        uint32_t hostTimeout;
        bool has_light;

        std::mutex connMutex;

        std::map< TcpConnection::Ptr, Host, PtrCompare<TcpConnection>> connections;
        std::map< TcpConnection::Ptr, std::set<std::string>, PtrCompare<TcpConnection>> subscriptions;
        
        boost::asio::io_service ioService;
        boost::thread_group pool;
        boost::asio::deadline_timer timer;
        std::shared_ptr<TcpServer> server;
        
        std::shared_ptr<Camera> cam;
        std::string outputdir;
        
        FrameQueue fqueue;
        Seeing seeing;
        Array<float> dd,ff,gg;
        float *ddPtr, *ggPtr;
        
        friend class TcpServer;

        
    };


}

#endif // WFWFS_DAEMON_HPP
