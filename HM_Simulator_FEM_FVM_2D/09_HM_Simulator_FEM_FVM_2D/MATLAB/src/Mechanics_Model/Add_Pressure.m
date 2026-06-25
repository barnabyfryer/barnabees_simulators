function [f] = add_pressure(f,Gen,P)

    for i = 1:Gen.Ne

        %Grab pressure of this element
        p = Gen.biot*P(i);

        %Find force to be applied to element edge, note only valid for
        %constnat pressure in element
        F_app_x = p * Gen.dy / 2;
        F_app_y = p * Gen.dx / 2;

        %Find nodes around element
        nodes = find(Gen.Ref(i,:) > 0);

        %Grab these nodes
        botleft_node = nodes(1);
        botright_node = nodes(2);
        topleft_node = nodes(3);
        topright_node = nodes(4);

        %Collect nodes on each side
        left_nodes = [botleft_node topleft_node];
        right_nodes = [botright_node topright_node];
        bot_nodes = [botleft_node botright_node];
        top_nodes = [topleft_node topright_node];

        %Apply appropriate force to each node group
        f(2*right_nodes-1) = f(2*right_nodes-1) + F_app_x;
        f(2*left_nodes-1) = f(2*left_nodes-1) - F_app_x;

        f(2*top_nodes) = f(2*top_nodes) + F_app_y;
        f(2*bot_nodes) = f(2*bot_nodes) - F_app_y;

    end
end
