#ifndef _champmanager_hpp
#define _champmanager_hpp

#include"champ.hpp"

class ChampManager
{

	private:
		std::vector<Champ*> lchamps_ ;
	public:
		ChampManager();
		~ChampManager();
		void initChamps(Config);
		void calculChamps(const Grid&,Sample&);
		void writeChamps(const Grid&);
		const Champ* getchamp(string) const; 
};
#endif
