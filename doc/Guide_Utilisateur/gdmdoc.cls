%=================================================
% gdmdoc.cls
% miniclasse to write very simple docs  
% Vincent Richefeu 2003, 2006
% last modifications : 15/02/2006
%.................................................
\ProvidesClass{simpledoc}
\LoadClass[twoside,a4paper]{refart}

\RequirePackage[french]{varioref}
\RequirePackage[french]{babel} 
%\RequirePackage[latin1]{inputenc} % pour PC
\RequirePackage[applemac]{inputenc} % pour MAC

\RequirePackage{amsmath}
\RequirePackage{graphicx}
\RequirePackage{float}
\RequirePackage[sort]{natbib}
\RequirePackage{tocloft}

%%
\RequirePackage{listings}

\lstdefinestyle{gdmcpp}{
language = C++,
frame = single,
numbers=left,
numberstyle=\tiny,
basicstyle=\small,
}

\lstdefinestyle{gdmdata}{
language = C++,
frame = single,
numbers=left,
numberstyle=\tiny,
basicstyle=\small,
emph = {Simulation,System,Algo,GroupData,GroupRelationData,Sample,Network} , emphstyle=\bfseries
}

\lstdefinestyle{mgpdata}{
language = C++,% a voir
frame = single,
numberstyle=\tiny,
basicstyle=\small
}

%%




