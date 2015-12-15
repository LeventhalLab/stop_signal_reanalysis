% script to keep track of the best wires to plot for each tetrode

IM164path = '/Volumes/dan/Recording data/IM-164_DL-22/IM-164_DL-22 screw turn hsds';
IM166path = '/Volumes/dan/Recording data/IM-166_DL-20/IM-166_DL-20 screw turn hsds';
turns_per_inch = 81;
turns_per_mm = turns_per_inch / 25.4;
mm_per_turns = 1/turns_per_mm;

IM164dates(1,:) = '2009-09-15';
IM164dates(2,:) = '2009-09-16';
IM164dates(3,:) = '2009-09-17';
IM164dates(4,:) = '2009-09-18';
IM164dates(5,:) = '2009-09-21';
IM164dates(6,:) = '2009-09-22';
IM164dates(7,:) = '2009-09-23';
IM164dates(8,:) = '2009-09-24';
IM164dates(9,:) = '2009-09-25';
IM164dates(10,:) = '2009-09-26';
IM164dates(11,:) = '2009-09-27';
IM164dates(12,:) = '2009-09-28';
IM164dates(13,:) = '2009-09-29';
IM164dates(14,:) = '2009-09-30';
IM164dates(15,:) = '2009-10-02';
IM164dates(16,:) = '2009-10-03';



IM166dates(1,:) = '2009-10-28';
IM166dates(2,:) = '2009-10-29';
IM166dates(3,:) = '2009-10-30';
IM166dates(4,:) = '2009-10-31';
IM166dates(5,:) = '2009-11-01';
IM166dates(6,:) = '2009-11-02';
IM166dates(7,:) = '2009-11-03';
IM166dates(8,:) = '2009-11-04';
IM166dates(9,:) = '2009-11-05';
IM166dates(10,:) = '2009-11-06';
IM166dates(11,:) = '2009-11-07';
IM166dates(12,:) = '2009-11-08';
IM166dates(13,:) = '2009-11-09';
IM166dates(14,:) = '2009-11-10';
IM166dates(15,:) = '2009-11-11';
IM166dates(16,:) = '2009-11-12';
IM166dates(17,:) = '2009-11-13';




%**************************************************************************
%**************************************************************************

i = 1;

tetrode{i}.IM = '164';
tetrode{i}.num = 1;

datevar = ['IM' tetrode{i}.IM 'dates'];
numDays = size(eval(datevar), 1);

for j = 1 : numDays
    
    evalstr = [datevar '(' num2str(j) ',:)'];
    tetrode{i}.dates(j,:) = eval(evalstr);
    
end

tetrode{i}.bestWire(1) = 1;
tetrode{i}.bestWire(2) = 1;
tetrode{i}.bestWire(3) = 1;
tetrode{i}.bestWire(4) = 1;
tetrode{i}.bestWire(5) = 1;
tetrode{i}.bestWire(6) = 1;
tetrode{i}.bestWire(7) = 1;
tetrode{i}.bestWire(8) = 1;
tetrode{i}.bestWire(9) = 1;
tetrode{i}.bestWire(10) = 1;
tetrode{i}.bestWire(11) = 1;
tetrode{i}.bestWire(12) = 1;
tetrode{i}.bestWire(13) = 1;
tetrode{i}.bestWire(14) = 1;
tetrode{i}.bestWire(15) = 1;
tetrode{i}.bestWire(16) = 1;
tetrode{i}.bestWire(17) = 1;

tetrode{i}.numTurns(1) = 10.0;
tetrode{i}.numTurns(2) = 2.0;
tetrode{i}.numTurns(3) = 0.25;
tetrode{i}.numTurns(4) = 1.25;
tetrode{i}.numTurns(5) = 0.5;
tetrode{i}.numTurns(6) = 0.125;
tetrode{i}.numTurns(7) = 0.125;
tetrode{i}.numTurns(8) = 0;
tetrode{i}.numTurns(9) = 0.0625;
tetrode{i}.numTurns(10) = 0.0625;
tetrode{i}.numTurns(11) = 0.0625;
tetrode{i}.numTurns(12) = 0.0625;
tetrode{i}.numTurns(13) = 0;
tetrode{i}.numTurns(14) = 0;
tetrode{i}.numTurns(15) = 0;
tetrode{i}.numTurns(16) = 0;

ML = 0;
DV = 0;
AP = 0;
tetrode{i}.pos(numDays, :) = [0 0 0];                % [ML DV AP]
for j = numDays - 1 : -1 : 1

    tetrode{i}.pos(j) = [ML, ...
                     tetrode{i}.pos(j+1,2) - (tetrode{i}.numTurns(j+1) * mm_per_turns), ...
                     AP];

end    % end for j...



