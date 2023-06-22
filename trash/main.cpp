// void process_node(dfg_node* n, int& curr_timestamp, int &tot_timestamps){
// 	// Put all of the next nodes in queue
// 	// Pop from queue,
// 	//	for each popped node:
// 	//		Check if node can be executed (i.e. lhs is done and rhs is done)
// 	//			if yes: execute queue, push next node in queue
// 	//			if no : push lhs or rhs or both into queue
	
// 	std::vector<dfg_node*> list;
// 	list.push_back(n);

// 	// n->process(curr_timestamp);

// 	dfg_node* curr_node;
// 	while(list.size() != 0)	{
// 		curr_node = list.front();

// 		for(int a = 0; a < list.size(); a++){
// 			std::cout << list[a]->label << " ";
// 		}
// 		std::cout << std::endl;
		
// 		// Check if node has already been executed
// 		if(curr_node->last_exec_time == curr_timestamp){
// 			pop_front(list);
// 			for(int i = 0; i < curr_node->next_nodes.size(); i++){
// 				list.push_back(curr_node->next_nodes[i]);
// 			}
// 		}
// 		else{
// 			// Check if node CAN be executed
// 			if(curr_node->t == CONST || curr_node->t == DELAY){
// 				std::cout << "[" << curr_timestamp << "] Processing " << curr_node->label << std::endl;
// 				curr_node->process(curr_timestamp);
// 				pop_front(list);
// 				for(int i = 0; i < curr_node->next_nodes.size(); i++){
// 					list.push_back(curr_node->next_nodes[i]);
// 				}
// 			}
// 			else if(curr_node->node_args_ready(curr_timestamp)){
// 				// Execute node
// 				std::cout << "[" << curr_timestamp << "] Processing " << curr_node->label << std::endl;
// 				curr_node->process(curr_timestamp);
// 				pop_front(list);
// 				for(int i = 0; i < curr_node->next_nodes.size(); i++){
// 					list.push_back(curr_node->next_nodes[i]);
// 				}
// 			}
// 			// Otherwise
// 			else{
// 				pop_front(list);
// 				// Schedule args execution
// 				if(curr_node->lhs != nullptr && curr_node->lhs->last_exec_time != curr_timestamp){
// 					// list.push_back(curr_node->lhs);
// 					list.insert(list.begin(), curr_node->lhs);
// 				}

// 				if(curr_node->rhs != nullptr && curr_node->rhs->last_exec_time != curr_timestamp){
// 					list.insert(list.begin(), curr_node->rhs);
// 					// list.push_back(curr_node->rhs);
// 				}

// 			}
// 		}
// 	}

// }


// void old_process_node(dfg_node* n, int& curr_timestamp, int &tot_timestamps);
// void old_process_node(dfg_node* n, int& curr_timestamp, int &tot_timestamps){
//     // If current node is processed, then process next node
    
//     // else, check if dependents are done processing,
//     //         if so, process current node

//     if(n->last_exec_time == curr_timestamp){
//         while(curr_timestamp < tot_timestamps){
//             for(int i = 0; i < n->next_nodes.size(); i++){
//                 if(n->next_nodes[i]->last_exec_time < curr_timestamp){
//                     process_node(n->next_nodes[i], curr_timestamp, tot_timestamps);
//                 }
//             }
//             curr_timestamp++;

//         }

//     }
//     else{
//         if( (n->t == DELAY && (curr_timestamp == n->last_exec_time + 1)) ||
// 			 n->t == CONST || 
// 			 n->t == INPUT )
// 		{
// 			std::cout << "[" << curr_timestamp << "] Processing " << n->label << std::endl;
//             n->process(curr_timestamp);
//         }
//         else{
//             for(int i = 0; i < n->prev_nodes.size(); i++){
//                 if(n->prev_nodes[i]->last_exec_time < curr_timestamp){
//                     process_node(n->prev_nodes[i], curr_timestamp, tot_timestamps);
//                 }
//             }
// 			std::cout << "[" << curr_timestamp << "] Processing " << n->label << std::endl;
//             n->process(curr_timestamp);
//         }
//     }

// }