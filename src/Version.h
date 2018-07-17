// Copyright (c) 2018 Michael C. Heiber & Joshua S. Brown
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCIMONTEC_VERSION_H
#define EXCIMONTEC_VERSION_H

#include <sstream>
#include <string>
#include <map>

namespace Excimontec {
	
  class Version {
    private:
      int major_ = 1;
      int minor_ = 0;
      int bug_release_ = 5;
      std::string release_type_ = "beta";
			static const std::map<std::string,int> key_type_ = {
				{"beta",0},
				{"alpha",1}};
		
    public:
			setVersion(std::string version){
				std::replace(version.begin(),version.end(),'v',' ');
				std::replace(version.begin(),version.end(),'.',' ');
				std::replace(version.begin(),version.end(),'-',' ');
				istringstream iss(version);
				iss >> major;
				iss >> minor;
				iss >> release_type_;
				iss >> bug_release_;
			}
      const int getMajor(){return major_;}
      const int getMinor(){return minor_;}
      const int getBugRelease(){return bug_release_;}
      const std::string getReleaseType(){return release_type_;}
      const std::string getVersion(){
        return "v"+to_string(major)+"."+to_string(minor)+"-"
          release_type_+"."+to_string(bug_release_);

			friend bool operator==(const Version &v1, const Version &v2);
			friend bool operator!=(const Version &v1, const Version &v2);
			friend bool operator>(const Version &v1, const Version &v2);
			friend bool operator>(const Version &v1, const Version &v2);
			friend bool operator<=(const Version &v1, const Version &v2);
			friend bool operator<=(const Version &v1, const Version &v2);

      }
  };

	// The default is the current version
	const Version current_version;

	bool operator==(const Version &v1, const Version &v2){
		if(v1.major_!=v2.major_) return false;
		if(v1.minor_!=v2.minor_) return false;
		if(v1.bug_release_!=v2.bug_release_) return false;
		if(v1.release_type_.compare(v2.release_type_)!=0) return false;
		return true;
	}

	bool operator!=(const Version &v1, const Version &v2){ return !(v1==v2);}
	
	bool operator<(const Version &v1, const Version &v2){
		if(v1.major_>v2.major_) return false;
		if(v1.minor_>v2.minor_) return false;
		if(key_type_[v1.release_type_]>key_type_[v2.release_type_]) return false;
		if(v1.bug_release_>v2.bug_release_) return false;
		if(v1==v2) return false;
		return true;
	}

	bool operator<=(const Version &v1, const Version &v2){
		if(v1==v2) return true;
		if(v1<v2) return true;
		return false;
	}

	bool operator>(const Version &v1, const Version &v2){
		if(v1<=v2) return false;
		return true;
	}

	bool operator>=(const Version &v1, const Version &v2){
		if(v1<v2) return true;
		return false;
	}
}

#endif // EXCIMONTEC_VERSION_H
