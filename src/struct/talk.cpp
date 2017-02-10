#include "talk.hpp"

void gdm::message(const char * txt)
{
  std::cout << " " << txt << std::endl << std::flush;
}

void gdm::warning(const char * txt)
{
  std::cout << " WARNING ! > " << txt << std::endl << std::flush;
}

void gdm::fatal(const char * txt)
{
  std::cerr << " FATAL ERROR ! > " << txt << std::endl << std::flush;
  exit(0);
}

void gdm::log(const char * txt, std::ostream & os)
{
  os << txt << std::endl << std::flush;
}

void gdm::fun_GDM(std::ostream & os)
{
  os << std::endl;
  os << "..................................................." << std::endl;
  os << " qjsdlakzi8   jksg&qsdusyè      d!çz       éjhsd   " << std::endl;
  os << "dskjfbiezué   dkjshcbemiumn     xépio     rjncxe   " << std::endl;
  os << "çà!             oir      yci    èyc oé   7èo ryj   " << std::endl;      
  os << "kqj             jks      xgq    sdu syè d!ç  éjh   " << std::endl;
  os << "dsk             jfb      iez    ué&  i!dkj   shc   " << std::endl;
  os << "çà!   urcoç     é5!      ryc    6oi   ryc    èyi   " << std::endl;
  os << "kqj   sdlak     zi8      73h    djk    T     sgg   " << std::endl;
  os << "dsk     j&b     iez      uéç    i!d          kjs   " << std::endl;
  os << "çà!     urc     oçé      è!r    yc8          oir   " << std::endl;
  os << "kqjsdlakzi8   o73hdjksggqsd    zusy          èdçmx " << std::endl;
  os << " dskjfbiez    uéçi!dkjshcb     eemi          umnwe " << std::endl;
  os << "................................... ToOL-KiT ......" << std::endl;
  os << std::endl << std::flush;
}

void gdm::debugMarker(const char* txt)
{
  std::cerr << "-> Marker: " << txt << std::endl << std::flush;
}


