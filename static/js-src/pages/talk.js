block post-footer
  script(type='html/template', id='relation-synopsis-template')
    p.lead.text-xs-center Click on the circles to see how your answers compared to the community's.
    ul#relation-synopsis-bar.list-unstyled.list-inline
    #feedback-next-action-area
    #chart-context(style='display:none;').row
      .col-xs-4.col-xs-offset-1.text-right
        p#concept-a.lead
      .col-xs-2.text-xs-center
        | <i class="fa fa-arrows-h fa-2x" aria-hidden="true"></i>
      .col-xs-4.text-left
        p#concept-b.lead
    ul#chart-list.list-unstyled

  script.
    var relation_task_settings = {
      'document_pk': "{{ doc.pk }}",
    };
  script(src="/static/js/tasks/relation-synopsis.js")

  script.
    var self_data, passages, regions;

    YPet.addInitializer(function(options) {

      $.getJSON('/task/entity-recognition/{{doc.pk}}/user/{{user.pk}}/results.json', function( data ) {
        self_data = data;
        passages = data.collection.document.passage;
        regions = {};

        _.each(passages, function(passage, passage_idx) {
          var passage_id = _.find(passage.infon, function(o){return o['@key']=='id';})['#text'];
          var p_body = '<div id="'+ passage_id +'" class="paragraph-box m-t-1"><p class="paragraph"></p></div></div>';
          $('.paragraphs').append(p_body);
          regions[''+passage_idx] = '#'+passage_id;
        });
        YPet.addRegions(regions);

        _.each(passages, function(passage, passage_idx) {
          var p = new Paragraph({'text': passage.text});
          YPet[''+passage_idx].show( new WordCollectionView({
            collection: p.get('words'),
            passage_json: passage,
            bioc_json: data
          }) );
          YPet[''+passage_idx].currentView.drawBioC(passage, false);
          YPet[''+passage_idx].currentView.drawBioC(null, true);
        });

      });
    });
    YPet.start();


    $('ul.pagination li a').on('mouseover', function(evt) {
      var user_pk = $(this).data('userpk');

      $.getJSON('/task/entity-recognition/{{doc.pk}}/user/'+ user_pk +'/results.json', function( data ) {
        self_data = data;
        passages = data.collection.document.passage;

        _.each(passages, function(passage, passage_idx) {
          var p = new Paragraph({'text': passage.text});
          YPet[''+passage_idx].show( new WordCollectionView({
            collection: p.get('words'),
            passage_json: passage,
            bioc_json: data
          }) );
          YPet[''+passage_idx].currentView.drawBioC(passage, false);
          YPet[''+passage_idx].currentView.drawBioC(null, true);
        });
      });

    });

