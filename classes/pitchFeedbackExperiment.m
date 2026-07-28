classdef pitchFeedbackExperiment
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        bird string
        date string
        numMotifs {isinteger}
        pitch (1,:) double
        isNoise (1,:) logical
        noiseLatency (1,:) double
        direction string
        noiseCorr double
        trigTimes (1,:) double
    end

    % methods
    %     function obj = untitled2(inputArg1,inputArg2)
    %         %UNTITLED2 Construct an instance of this class
    %         %   Detailed explanation goes here
    %         obj.Property1 = inputArg1 + inputArg2;
    %     end
    % 
    %     function outputArg = method1(obj,inputArg)
    %         %METHOD1 Summary of this method goes here
    %         %   Detailed explanation goes here
    %         outputArg = obj.Property1 + inputArg;
    %     end
    % end
end