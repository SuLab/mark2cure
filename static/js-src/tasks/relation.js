
$('#next_button').on('click', function(evt) {
  /* Causes the next to load */
  var current_relationship = collection.findWhere({'current': true});
  var next_relationship = collection.findWhere({user_completed: false, available: true});

  if (next_relationship) {
    current_relationship.set('user_completed', true);
    $('#feedback-next-action-area').hide();
    $('#tree-action-area').show();
    $('#chart').empty();
    $('#chart-list').empty();
  } else {
    $('#tree-action-area').hide();
    $('#task_relation_submit_document_set').submit();
  }

});

