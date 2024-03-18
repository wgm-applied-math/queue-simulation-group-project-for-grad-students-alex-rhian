classdef Renege < Event
    % Departure Subclass of Event that represents the departure of a
    % Customer.

    properties

        % ID of the customer that reneges- their alarm has gone off!
        % remove from queue if they have NOT been served, if they are being
        % served or have already been served then do nothing 

        % ServerIndex - Index of the service station from which the
        % departure occurred
        % We need to change this to check the person's idea to see if when
        % their alarm goes off, are they still in line?

        Id;
    end
    methods
        function obj = Renege(Time, Id)

            % Departure - Construct a renege event from a time and
            % customer Id.

            % Departure - Construct a departure event from a time and
            % server index.
 
            arguments
                Time = 0.0;
                Id = 0;
            end
            
            % MATLAB-ism: This incantation is how to invoke the superclass
            % constructor.
            obj = obj@Event(Time);

            obj.Id = Id;
        end
        function varargout = visit(obj, other)
            % visit - Call handle_departure

            % MATLAB-ism: This incantation means whatever is returned by
            % the call to handle_departure is returned by this visit
            % method.
            [varargout{1:nargout}] = handle_renege(other, obj);
        end
    end
end