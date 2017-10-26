var channel = Backbone.Radio.channel('mark2cure');

//
// Models + Collections
//

TalkDocumentNERAnn = Backbone.Model.extend({
  defaults: {
    'text': '',
    'occurances': 0
  }
});

TalkDocumentContributor = Backbone.Model.extend({
  defaults: {
    user_id: null
  }
});

TalkDocumentContributorList = Backbone.Collection.extend({
  model: TalkDocumentContributor,
  url: function() {
    return '/api/talk/document/'+this.document_pk+'/'+this.task_type+'/contributor/list/'
  },
  initialize: function(options) {
    this.document_pk = options['document_pk']
    this.task_type = 'ner'
  }
})

TalkDocumentNERAnnList = Backbone.Collection.extend({
  model: TalkDocumentNERAnn,
  url: function() {
    return '/api/talk/document/'+this.document_pk+'/annotations/'+this.ann_type_idx+'/list/'
  },

  initialize: function(options) {
    this.document_pk = options['document_pk'];
    this.ann_type_idx = options['ann_type_idx'];
  }
})

TalkComment = Backbone.Model.extend({
  defaults: {
    'user_id': null,
    'user_name': '',
    'comment': '',
    'submit_date': null,
    'document_pk': null
  }
});


TalkCommentList = Backbone.Collection.extend({
  model: TalkComment,
  url: '/api/talk/comment/list/'
});


TalkDocument = Backbone.Model.extend({
  defaults: {
    'id': null,
    'title': '',
    'comments': 0
  }
});


TalkDocumentList = Backbone.Collection.extend({
  model: TalkDocument,
  url: '/api/talk/document/list/'
});

//
// Views
//
TalkDocumentItemView = Backbone.Marionette.View.extend({
  template: '#talk-document-item-template',
  className: 'list-group-item flex-column align-items-start',
});

TalkDocumentListView = Backbone.Marionette.CollectionView.extend({
  childView: TalkDocumentItemView
});

TalkDiscussedDocumentView = Backbone.Marionette.View.extend({
  template: '#talk-list-template',
  templateContext: {
    'list_name': 'Discussed Documents'
  },
  className: 'col',

  regions: {
    'list': '.list-group',
  },

  collectionEvents: {
    'sync': function() { this.render(); }
  },

  initialize: function() {
    this.collection = new TalkDocumentList();
    this.collection.fetch();
  },

  onRender: function() {
    this.showChildView('list', new TalkDocumentListView({'collection': this.collection}));
  }
});

TalkCommentItemView = Backbone.Marionette.View.extend({
  template: '#talk-comment-item-template',
  className: 'list-group-item flex-column align-items-start',
  ui: {
    'time': '.time-ago',
  },
  onRender: function() {
    this.ui.time.html( moment.utc(this.model.get('submit_date')).format('LL') );
  }
});

TalkCommentListView = Backbone.Marionette.CollectionView.extend({
  childView: TalkCommentItemView
});

TalkCommentView = Backbone.Marionette.View.extend({
  template: '#talk-list-template',
  templateContext: {
    'list_name': 'Recent Comments'
  },
  className: 'col',

  regions: {
    'list': '.list-group',
  },

  modelEvents: {
    'sync': function() { this.render(); }
  },

  initialize: function() {
    this.collection = new TalkCommentList();
    this.collection.fetch();
  },

  onRender: function() {
    this.showChildView('list', new TalkCommentListView({'collection': this.collection}));
  }
});


TalkView = Backbone.Marionette.View.extend({
  template: '#talk-template',
  className: 'row justify-content-center mt-3',

  regions: {
    'documents': '#talk-documents',
    'comments': '#talk-comments'
  },

  onRender: function() {
    this.showChildView('documents', new TalkDiscussedDocumentView());
    this.showChildView('comments', new TalkCommentView());
  }
});

TalkDocumentNERUserItem = Backbone.Marionette.View.extend({
  template: '#talk-document-user-item-template',
  tagName: 'li',
  className: 'list-inline-item',
  events: {
    'mouseover': function(evt) {
      $.getJSON('/task/ner/'+ channel.request('get:document:pk') +'/user/'+ this.model.get('user_id')  +'/results.json', function( data ) {
        console.log(data);
      });
    }
  }
});

TalkDocumentNERUserList = Backbone.Marionette.CollectionView.extend({
  childView: TalkDocumentNERUserItem,
  tagName: 'ul',
  className: 'list-inline'
});

TalkDocumentNERView = Backbone.Marionette.View.extend({
  template: '#talk-document-ner-template',
  className: 'row justify-content-center',

  regions: {
    'contributors': '#talk-ner-user-list',
    'ner_viz': '#talk-ner-display'
  },

  collectionEvents: {
    'sync': function() { this.render(); }
  },

  initialize: function() {
    this.collection = new TalkDocumentContributorList(this.options)
    this.collection.fetch();
  },

  onRender: function() {
    this.showChildView('contributors', new TalkDocumentNERUserList({'collection': this.collection}));
  }
});


TalkDocumentView = Backbone.Marionette.View.extend({
  template: '#talk-document-template',
  className: 'row justify-content-center mt-3',

  regions: {
    'ner': '#talk-document-ner',
    'comment': '#talk-document-comment',
  },

  initialize: function() {
    channel.reply('get:document:pk', this.getDocumentPk, this);
  },

  onRender: function() {
    // this.collection = new TalkDocumentNERAnnList({'document_pk':4229, 'ann_type_idx':1});
    // this.collection.fetch();
    this.showChildView('ner', new TalkDocumentNERView(this.options));
  },

  getDocumentPk: function() {
    return this.getOption('document_pk');
  }
});

// block post-footer
//   script(type='html/template', id='relation-synopsis-template')
//     p.lead.text-xs-center Click on the circles to see how your answers compared to the community's.
//     ul#relation-synopsis-bar.list-unstyled.list-inline
//     #feedback-next-action-area
//     #chart-context(style='display:none;').row
//       .col-xs-4.col-xs-offset-1.text-right
//         p#concept-a.lead
//       .col-xs-2.text-xs-center
//         | <i class="fa fa-arrows-h fa-2x" aria-hidden="true"></i>
//       .col-xs-4.text-left
//         p#concept-b.lead
//     ul#chart-list.list-unstyled
//
//   script.
//     var relation_task_settings = {
//       'document_pk': "{{ doc.pk }}",
//     };
//   script(src="/static/js/tasks/relation-synopsis.js")
//
//   script.
//     var self_data, passages, regions;
//
//     YPet.addInitializer(function(options) {
//
//       $.getJSON('/task/entity-recognition/{{doc.pk}}/user/{{user.pk}}/results.json', function( data ) {
//       $.getJSON('/task/entity-recognition/{{doc.pk}}/user/'+ user_pk +'/results.json', function( data ) {
