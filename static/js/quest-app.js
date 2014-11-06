var activate_tabs = function(doc_id) {
  $('.guest-guide .doc-item-col').removeClass('active');
  var $tab_el = $('.quest-guide').find("[data-doc='" + doc_id + "']");
  $tab_el.addClass('active');
  $tab_el.prevAll().addClass('muted').removeClass('active');

  $('#insert .container.document').hide();
  $('#doc_'+doc_id).fadeIn();
};

$('#quest-submit').on('click', function(evt) {
  /* Options:
   * 1. Submit Doc
   * 2. Go to next Doc
   * 3. Complete / Submit Quest
   */
  var $document = $("#insert .container.document:visible");
  var document_id = $document.data('doc');
  var $document_form = $document.find('form');
  var sections = $document.find('.game div').map(function(k,v) { return +$(v).attr('id'); });

  $(this).attr('disabled', 'disabled');
  /* Submit the annotations to the server, then tell the server
  you're done with that document. */

  /* We prevent this or else Django gives broken pipe b/c we need to
  wait for the ajax submissions to close before reloading the page */
  evt.preventDefault();

  if($document.find('.results').length) {
    var doc_ids = [];
    $('.doc-item-col').each(function() { doc_ids.push($(this).data('doc')); });

    if( doc_ids.indexOf(document_id) == doc_ids.length-1) {
      /* Submit the Quest */
      $('#quest-complete').submit();
    } else {
      /* Show the next Document */
      var next_id = doc_ids[ doc_ids.indexOf(document_id)+1 ];
      activate_tabs(next_id);
    };

    /* Currently looking at results */
    $(this).text('Submit');
    $(this).attr('disabled', false);
  } else {
    var counter = 0,
        ann_counter = 0;

    /* Iterate over each of the paragraphs or annotatable sections on the page */
    _.each(sections, function(section_id) {
      var annotations = YPet[section_id].currentView.collection.parentDocument.get('annotations').toJSON(),
          url = '/document/'+ task_id +'/'+ document_id +'/section/'+ section_id +'/annotation/create/';
      ann_counter += annotations.length;

      var csrf = $document.find('input[name="csrfmiddlewaretoken"]').val();

      /* Iterate over each of the annotations within that section */
      _.each(annotations, function(annotation) {
        $.ajax({
          type: 'POST',
          url: url,
          data: $.extend({'csrfmiddlewaretoken': csrf}, annotation),
          cache: false,
          async: false,
          success: function() { counter++; },
        });
      });
    });

    /* If they all got sent to the server, let's move on */
    if(counter === ann_counter) {
      var form = $document_form,
          data = {},
          input;

      $.each( form.find('input, textarea, select'), function() {
        input = $(this);
        data[ input.attr('name') ] = input.val();
        input.val('')
      });

      $.ajax({
        type: 'POST',
        url: form.attr('action'),
        dataType: 'html',
        data: data,
        cache: false,
        async: false,
        success: function(data) {
          $document.html(data);
          $pag_tab = $('.pagination li[data-doc='+document_id+']');
          $pag_tab.addClass('disabled');
          $pag_tab.next().addClass('active').trigger('click');
        }
      });

      $(this).text('Next');
      $(this).attr('disabled', false);
    } else {
      alert('There was a problem submissing this document. Please try another.');
    }
  };

});
