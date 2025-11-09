function blkStruct = slblocks
%SLBLOCKS Defines the block library for NCD Blockset.

%   $Revision: 1.3 $ $Date: 1998/08/07 17:51:25 $
%   Copyright (c) 1990-1998 by The MathWorks, Inc. All Rights Reserved.

% Name of the subsystem which will show up in the
% SIMULINK Blocksets and Toolboxes subsystem.
blkStruct.Name = 'lucas';

% The function that will be called when
% the user double-clicks on this icon.
blkStruct.OpenFcn = 'lucas';

% The argument to be set as the Mask Display for the subsystem.
% You may comment this line out if no specific mask is desired.
blkStruct.MaskDisplay = ['plot([0:10],[-.5 1.5 .6 1.3 .8 1.1 .95 1.02 .99 1 1]);disp(''\n NCD'')'];

% Define the library list for the Simulink Library browser.
% Return the name of the library model and the name for it
Browser(1).Library = 'lucas';
Browser(1).Name    = 'lucas';

blkStruct.Browser = Browser;


% [End] of slblocks.m
