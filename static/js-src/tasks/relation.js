var submit_status = function() {
  /* Submit button status changes accordingly */
  if( RelationApp.start.currentView.options.choice ) {
    $('#submit_button').attr('disabled', false).removeClass('disabled');
  } else {
    $('#submit_button').attr('disabled', true).addClass('disabled');
  };
};




$('#next_button').on('click', function(evt) {
  /* Causes the next to load */
  var current_relationship = collection.findWhere({'current': true});
  var next_relationship = collection.findWhere({user_completed: false, available: true});

  if (next_relationship) {
    current_relationship.set('user_completed', true);
    submit_status();
    $('#feedback-next-action-area').hide();
    $('#tree-action-area').show();
    $('#chart').empty();
    $('#chart-list').empty();
  } else {
    $('#tree-action-area').hide();
    $('#task_relation_submit_document_set').submit();
  }

});

$('#submit_button').on('click', function(evt) {
  /* When a user is confirming their Tree choice selection */

  var current_selection = RelationApp.start.currentView.options.choice;
  var current_relationship = collection.findWhere({'current': true});

  if(current_selection.get('id')) {
    /* Submit the data from the current selection to the server */
    $.ajax({
      type: 'POST',
      url: '/task/relation/'+ relation_task_settings.document_pk +'/'+ current_relationship.id +'/submit/',
      data: $.extend({'csrfmiddlewaretoken': relation_task_settings.csrf_token },
                     {'relation': current_selection.get('id')}),
      cache: false,
      success: function() {

        /* If they flag a concept as incorrect, remove it from the
         * availability pool of tasks for this relationship document */
        if( _.contains(incorrect_flag_arr, current_selection.get('id')) ) {
          var incorrect_key = 'c'+(1+incorrect_flag_arr.indexOf(current_selection.get('id')));
          var incorrect_concept = current_relationship.get('concepts')[ incorrect_key ];
          collection.each(function(relation) {
            var concepts = relation.get('concepts');
            if( concepts['c1'].id == incorrect_concept['id'] || concepts['c2'].id == incorrect_concept['id'] ) {
              relation.set('available', false);
            }
          });
        }

        /* Have them review the community's answers */
        show_results(current_relationship.get('document_id'), current_relationship.get('id'));
        update_score();

      },
      error: function() { alert('Please refresh your browser and try again.') },
    });
  }
});
